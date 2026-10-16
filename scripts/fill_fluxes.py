#!/usr/bin/env python
"""
Gap-fill momentum and turbulent (latent/sensible) heat fluxes in a GLEN
quality-controlled half-hourly NetCDF file.

Methods
-------
friction_velocity : observed -> time-interp (<=3h) -> COARE u* solved from wind
momentum_flux     : observed -> rho*u*^2 (rho = P/(Rd*T*(1+0.608 q))) -> interp(<=3h)
latent/sensible   : observed -> time-interp (<=3h) -> Monin-Obukhov bulk
                    (Charusombat 2018) wherever the mooring Tw is available

Original observed values are never altered. A companion <var>_fill_method int8
flag variable records the provenance of every point.  Gaps with no surface
temperature for the MO calculation are left as NaN.

Usage
-----
    python fill_fluxes.py [input.nc] [output.nc]

Defaults to the Stannard Rock file.  Reusable for the other GLEN stations.
"""
import sys
import os
import numpy as np
import pandas as pd
import netCDF4 as nc
from datetime import datetime, timezone
from mo_flux import compute_mo_fluxes, compute_mo_fluxes_wind, ustar_from_wind
from mooring_lst import lst_at

# ---------------------------------------------------------------- constants
RD = 287.05          # J kg-1 K-1, gas constant for dry air
EPS = 0.608          # virtual-temperature coefficient
SHORT_GAP_H = 3.0    # max gap (hours) eligible for linear time interpolation
USTAR_MIN = 0.0      # physical floor for friction velocity
# PCHIP-fill gap limits on the half-hourly grid (48 samples/day)
SAMPLES_PER_DAY = 48
PRES_PCHIP_STEPS = SAMPLES_PER_DAY * 8   # air_pressure: interior gaps <= 8 days
WIND_PCHIP_STEPS = SAMPLES_PER_DAY * 1   # wind vector: interior gaps <= 24 h

# Monin-Obukhov winter heat-flux fill (Charusombat et al. 2018, via mo_flux.py)
# Lake-surface temperature comes from the Eastern Mooring near-surface sensor.
# MO_METHOD selects the driver for the bulk turbulent heat fluxes:
#   "wind"  -> full bulk algorithm driven by the wind speed at the measurement
#              height (u* solved internally), exactly as in Charusombat 2018
#              (default; verified to match the eddy-covariance fluxes better).
#   "ustar" -> flux-profile form using the measured/filled friction velocity u*
# Override at runtime with the GLEN_MO_METHOD environment variable.
MO_METHOD = os.environ.get("GLEN_MO_METHOD", "wind").lower()
# MO_COOL_SKIN applies the COARE cool-skin correction (Fairall et al. 1996) so
# the MO fluxes use the interface temperature (Tw - dter) rather than the bulk
# water temperature.  It needs downwelling shortwave & longwave; it is silently
# disabled at stations/timesteps where those are unavailable.  Override with the
# GLEN_MO_COOL_SKIN environment variable ("1"/"0").
MO_COOL_SKIN = os.environ.get("GLEN_MO_COOL_SKIN", "1") not in ("0", "false", "False")
WINDOW_DAYS = 60     # length of each winter window from its onset date
WINTER_START_CSV = (
    "/eagle/projects/COMPASS-GLM/wpringle/Lake_Julia_Runs/"
    "V1_lake_column_model_LES_vs_TURB/NORiOceanParameterization.jl/"
    "figure_data/lake_superior_eastern_mooring/"
    "lake_superior_eastern_mooring_winter_start_dates.csv")

# fill-method flag codes (shared vocabulary; per-variable subset applies).
# Code 3 (former random-forest fill) is retired and intentionally skipped so
# the remaining codes keep their meaning vs. earlier outputs.
(FLAG_OBS, FLAG_INTERP, FLAG_RHOU2, FLAG_WINDBIN, FLAG_MO, FLAG_COARE,
 FLAG_CLIM, FLAG_MAGNUS, FLAG_PCHIP, FLAG_VECTOR) = 0, 1, 2, 4, 5, 6, 7, 8, 9, 10
FLAG_MISSING = -1
FLAG_FILLVALUE = np.int8(-1)
FLAG_VALUES = np.array([0, 1, 2, 4, 5, 6, 7, 8, 9, 10], dtype=np.int8)
FLAG_MEANINGS = ("observed time_interpolated computed_rho_ustar_squared "
                 "wind_speed_bin_median monin_obukhov_bulk "
                 "coare_ustar_from_wind monthly_climatology "
                 "specific_humidity_magnus pchip_interpolated "
                 "wind_vector_interpolated")
_MO_DRIVER_DESC = ("measured/filled u*" if MO_METHOD == "ustar"
                   else "wind speed at measurement height (u* solved)")
FLAG_COMMENT = ("0=observed; 1=time-interpolated (<=%gh gap); "
                "2=computed rho*u*^2; 4=wind-speed-bin/overall median; "
                "5=Monin-Obukhov bulk (%s, Eastern-Mooring near-surface Tw); "
                "6=COARE u* solved from wind speed (Charusombat 2018; "
                "stability where winter Tw known, else neutral); "
                "7=monthly-mean climatology (air_pressure long gaps); "
                "8=specific_humidity recomputed from T,RH,P (Magnus); "
                "9=PCHIP interpolation (air_pressure gaps <=8 days); "
                "10=wind vector (u,v) PCHIP interpolation "
                "(wind_speed/direction, <=24h gap, cyclical); "
                "-1=still missing" % (SHORT_GAP_H, _MO_DRIVER_DESC))


# ---------------------------------------------------------------- helpers
def getvar(ds, name):
    """Return variable data as float64 ndarray with masked/fill -> NaN."""
    v = ds.variables[name][:]
    a = np.ma.filled(np.ma.masked_array(v, dtype="float64"), np.nan).astype("float64")
    return np.asarray(a)


def getvar_opt(ds, name):
    """Like getvar but returns None if the variable is absent (some stations
    lack downwelling shortwave/longwave)."""
    return getvar(ds, name) if name in ds.variables else None


def interp_short_gaps(t, y, max_gap_h):
    """Linear-interpolate NaNs in y vs time t, but only across gaps whose
    total elapsed time <= max_gap_h hours. Returns (filled, newly_filled_mask).
    Interior points only (no extrapolation beyond first/last valid sample)."""
    y = y.copy()
    valid = np.isfinite(y)
    if valid.sum() < 2:
        return y, np.zeros_like(y, dtype=bool)
    # full linear interpolation over interior region (np.interp clamps ends)
    yi = np.interp(t, t[valid], y[valid], left=np.nan, right=np.nan)
    # measure the elapsed-time span of each contiguous NaN run and the
    # bracketing valid times so we only accept short gaps
    fill_mask = np.zeros_like(y, dtype=bool)
    miss = ~valid
    idx = np.where(miss)[0]
    if idx.size == 0:
        return y, fill_mask
    runs = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
    for run in runs:
        i0, i1 = run[0], run[-1]
        if i0 == 0 or i1 == len(y) - 1:
            continue  # edge gap -> no interior bracket -> skip (no extrapolation)
        span_h = (t[i1 + 1] - t[i0 - 1]) / 3600.0
        if 0 < span_h <= max_gap_h:
            sel = run[np.isfinite(yi[run])]
            y[sel] = yi[sel]
            fill_mask[sel] = True
    return y, fill_mask


def specific_humidity_magnus(T_K, RH_pct, P_Pa):
    """Specific humidity (kg/kg) from temperature, relative humidity, pressure
    using the Magnus saturation-vapour-pressure relation over water -- the same
    method noted in the file for the original specific_humidity variable."""
    Tc = T_K - 273.15
    es = 6.112 * np.exp(17.62 * Tc / (243.12 + Tc))      # hPa
    e = (RH_pct / 100.0) * es                            # hPa
    Phpa = P_Pa / 100.0
    return 0.622 * e / (Phpa - 0.378 * e)                # kg/kg


def winter_window_mask(t_seconds):
    """Boolean mask of timesteps within WINDOW_DAYS of each winter-onset date
    (read from WINTER_START_CSV).  Used to scope the Monin-Obukhov heat-flux
    fill to the open-water mixing period where Tw~4C is defensible.  Returns an
    all-False mask if the CSV is unavailable."""
    mask = np.zeros(t_seconds.shape, dtype=bool)
    if not os.path.exists(WINTER_START_CSV):
        print(f"  WARNING: winter-start CSV not found, MO fill disabled:\n    "
              f"{WINTER_START_CSV}")
        return mask
    span = WINDOW_DAYS * 86400.0
    with open(WINTER_START_CSV) as f:
        header = f.readline().strip().split(",")
        si = header.index("start_date")
        for line in f:
            parts = line.strip().split(",")
            if len(parts) <= si or not parts[si]:
                continue
            d = datetime.strptime(parts[si], "%Y-%m-%d").replace(tzinfo=timezone.utc)
            s0 = d.timestamp()
            mask |= (t_seconds >= s0) & (t_seconds < s0 + span)
    return mask


def month_index(t_seconds):
    """Calendar month (1..12) for each epoch-second timestamp (UTC)."""
    out = np.empty(t_seconds.shape, dtype=int)
    for i, ts in enumerate(t_seconds):
        if np.isfinite(ts):
            out[i] = datetime.fromtimestamp(float(ts), tz=timezone.utc).month
        else:
            out[i] = -1
    return out


# ---------------------------------------------------------------- main
def main():
    infile = sys.argv[1] if len(sys.argv) > 1 else \
        "US_StannardRockSuperior_processed_halfhourly_qc.nc"
    outfile = sys.argv[2] if len(sys.argv) > 2 else \
        infile.replace(".nc", "_gapfilled.nc")

    print(f"Input : {infile}")
    print(f"Output: {outfile}\n")

    ds = nc.Dataset(infile)
    t = getvar(ds, "time")

    # ---- load raw target + predictor fields ------------------------------
    ust_raw = getvar(ds, "friction_velocity")
    mfl_raw = getvar(ds, "momentum_flux")
    le_raw = getvar(ds, "latent_heat_flux")
    sh_raw = getvar(ds, "sensible_heat_flux")

    T = getvar(ds, "air_temperature")
    U = getvar(ds, "wind_speed")
    WD = getvar(ds, "wind_direction")
    q = getvar(ds, "specific_humidity")
    RH = getvar(ds, "relative_humidity")
    P = getvar(ds, "air_pressure")
    Rs = getvar_opt(ds, "downwelling_shortwave_flux")   # for cool-skin (optional)
    Rl = getvar_opt(ds, "downwelling_longwave_flux")

    n = t.size
    print(f"timesteps: {n}")

    # =====================================================================
    # METEOROLOGICAL DRIVERS: short-gap interpolation (+ special long-gap fills
    # for pressure and humidity).  These filled series both drive the flux
    # calculations AND are written to the output with their own *_fill_method
    # flags.  Observed values are never altered; only gaps are filled.
    # =====================================================================
    def fill_interp(y):
        """Short-gap interpolation; returns (filled, int8 flag)."""
        flag = np.where(np.isfinite(y), FLAG_OBS, FLAG_MISSING).astype(np.int8)
        yf, m = interp_short_gaps(t, y, SHORT_GAP_H)
        flag[m] = FLAG_INTERP
        return yf, flag

    Tf, T_flag = fill_interp(T)
    qf, q_flag = fill_interp(q)
    RHf, RH_flag = fill_interp(RH)

    # wind_speed + wind_direction: interpolate the wind VECTOR (u,v components)
    # so the cyclical direction wraps correctly (rather than linearly
    # interpolating direction across the 0/360 discontinuity).  PCHIP fills
    # interior gaps up to 24 h (limit_area="inside" prevents edge extrapolation).
    # Speed and direction are recovered from the interpolated components;
    # observed values are preserved exactly. (wind_from_direction convention.)
    th = np.deg2rad(WD)
    u_wind = -U * np.sin(th)
    v_wind = -U * np.cos(th)
    uf = pd.Series(u_wind).interpolate(method="pchip", limit=WIND_PCHIP_STEPS,
                                       limit_area="inside").values
    vf = pd.Series(v_wind).interpolate(method="pchip", limit=WIND_PCHIP_STEPS,
                                       limit_area="inside").values
    spd = np.hypot(uf, vf)
    drc = np.degrees(np.arctan2(-uf, -vf)) % 360.0
    Uf = U.copy()
    WDf = WD.copy()
    gap_U = ~np.isfinite(U)
    gap_WD = ~np.isfinite(WD)
    Uf[gap_U] = spd[gap_U]
    WDf[gap_WD] = drc[gap_WD]
    U_flag = np.where(np.isfinite(U), FLAG_OBS, FLAG_MISSING).astype(np.int8)
    WD_flag = np.where(np.isfinite(WD), FLAG_OBS, FLAG_MISSING).astype(np.int8)
    U_flag[gap_U & np.isfinite(Uf)] = FLAG_VECTOR
    WD_flag[gap_WD & np.isfinite(WDf)] = FLAG_VECTOR

    # air_pressure: PCHIP interpolation across interior gaps up to ~8 days
    # (half-hourly data -> 48*8 steps), then monthly-mean climatology for longer
    # gaps and the leading/trailing edges.  limit_area="inside" keeps PCHIP from
    # extrapolating past the first/last observation (which blows up badly).
    P_flag = np.where(np.isfinite(P), FLAG_OBS, FLAG_MISSING).astype(np.int8)
    Pf = pd.Series(P).interpolate(method="pchip", limit=PRES_PCHIP_STEPS,
                                  limit_area="inside").values
    m_pchip = (~np.isfinite(P)) & np.isfinite(Pf)
    P_flag[m_pchip] = FLAG_PCHIP

    # downwelling radiation for the cool-skin correction (short-gap interp only;
    # remaining gaps are left NaN -> cool-skin auto-disables at those points)
    cool_ok = MO_COOL_SKIN and (Rs is not None) and (Rl is not None)
    if cool_ok:
        Rsf, _ = interp_short_gaps(t, Rs, SHORT_GAP_H)
        Rlf, _ = interp_short_gaps(t, Rl, SHORT_GAP_H)
    else:
        Rsf = Rlf = None
        if MO_COOL_SKIN and (Rs is None or Rl is None):
            print("  cool-skin requested but downwelling SW/LW not both "
                  "available; cool-skin disabled for this station")

    # remaining pressure gaps -> monthly-mean climatology (from OBSERVED
    # pressure, so PCHIP fills can't bias it), then overall observed mean
    months = month_index(t)
    if np.isnan(Pf).any():
        overall_P = np.nanmean(P)
        for m in range(1, 13):
            sel = (months == m) & np.isnan(Pf)
            if sel.any():
                clim = np.nanmean(P[months == m])
                Pf[sel] = clim if np.isfinite(clim) else overall_P
                P_flag[sel] = FLAG_CLIM
        rem = np.isnan(Pf)
        Pf[rem] = overall_P
        P_flag[rem] = FLAG_CLIM
    print(f"pressure filled for rho: {np.isfinite(Pf).sum()}/{n}")

    # recompute specific humidity (Magnus) where missing but T & RH available.
    # q was originally derived from T, RH, P; most q gaps coincide with P gaps,
    # so reconstructing it (with filled P) unlocks the flux fills that need it.
    need_q = (~np.isfinite(qf)) & np.isfinite(Tf) & np.isfinite(RHf) & np.isfinite(Pf)
    # sanity: agreement with original q where both exist
    chk = np.isfinite(q) & np.isfinite(Tf) & np.isfinite(RHf) & np.isfinite(Pf)
    if chk.any():
        qest = specific_humidity_magnus(Tf[chk], RHf[chk], Pf[chk])
        bias = np.nanmean(qest - q[chk]); rms = np.sqrt(np.nanmean((qest - q[chk])**2))
        print(f"q Magnus check vs original: bias={bias:.2e} rms={rms:.2e} kg/kg (n={chk.sum()})")
    qf[need_q] = specific_humidity_magnus(Tf[need_q], RHf[need_q], Pf[need_q])
    q_flag[need_q] = FLAG_MAGNUS
    print(f"specific_humidity recomputed at {need_q.sum()} points; "
          f"now available {np.isfinite(qf).sum()}/{n}")

    # Eastern-Mooring near-surface water temperature + winter windows; used both
    # for the COARE u* stability (section 1) and the Monin-Obukhov heat-flux
    # fill (sections 3-4).  Tw_winter restricts the mooring Tw to the well-mixed
    # winter windows where it is a defensible surface temperature.
    winter_mask = winter_window_mask(t)
    Tw_arr = lst_at(t)                       # degC at flux times, NaN if no data
    Tw_winter = np.where(winter_mask & np.isfinite(Tw_arr), Tw_arr, np.nan)

    # Southern-Mooring near-surface water temperature, for a sensitivity check
    # on the stand-alone COARE-from-wind diagnostics (section below) only --
    # the primary gap-fill above continues to use the Eastern-Mooring Tw.
    Tw_arr_south = lst_at(t, station="southern")
    n_tw_south = int(np.isfinite(Tw_arr_south).sum())
    print(f"Southern-Mooring Tw available at {n_tw_south}/{n} timesteps "
          f"(for COARE-from-wind sensitivity diagnostic)")

    # =====================================================================
    # 1. FRICTION VELOCITY
    # =====================================================================
    ust = ust_raw.copy()
    ust_flag = np.where(np.isfinite(ust), FLAG_OBS, FLAG_MISSING).astype(np.int8)

    # 1a. short-gap interpolation
    ust, m_interp = interp_short_gaps(t, ust, SHORT_GAP_H)
    ust_flag[m_interp] = FLAG_INTERP

    # 1b. COARE u* solved from wind speed (Charusombat 2018) for the remainder.
    #     Stability uses the winter-window mooring Tw where known, else neutral.
    ust_coare, _ = ustar_from_wind(Uf, Tf, qf, Pf, Tw_C=Tw_winter)
    use_coare = (~np.isfinite(ust)) & np.isfinite(ust_coare)
    ust[use_coare] = ust_coare[use_coare]
    ust_flag[use_coare] = FLAG_COARE

    # 1c. residual (e.g. wind missing) -> overall observed-u* median fallback
    obs = np.isfinite(ust_raw)
    if obs.sum() > 0:
        overall_med_ust = np.median(ust_raw[obs])
        resid = ~np.isfinite(ust)
        ust[resid] = overall_med_ust
        ust_flag[resid] = FLAG_WINDBIN
    ust = np.clip(ust, USTAR_MIN, None)
    print(f"u* filled: COARE-from-wind {int(use_coare.sum())}, "
          f"median fallback {int((ust_flag == FLAG_WINDBIN).sum())}")

    # =====================================================================
    # density (moist air, virtual temperature)
    # =====================================================================
    # q affects rho only through the ~1% virtual-temperature term; where q is
    # still unavailable but T is, fall back to dry air (q=0) -- the user noted
    # mom_flux = dry_air * u*^2 gives equivalent values. rho is then defined
    # wherever T (and the always-filled P) is available.
    q_for_rho = np.where(np.isfinite(qf), qf, 0.0)
    Tv = Tf * (1.0 + EPS * q_for_rho)
    rho = Pf / (RD * Tv)   # kg m-3 ; finite wherever Tf is finite

    # =====================================================================
    # 2. MOMENTUM FLUX  tau = rho * u*^2
    # =====================================================================
    mfl = mfl_raw.copy()
    mfl_flag = np.where(np.isfinite(mfl), FLAG_OBS, FLAG_MISSING).astype(np.int8)

    can_compute = (~np.isfinite(mfl)) & np.isfinite(ust) & np.isfinite(rho)
    mfl[can_compute] = rho[can_compute] * ust[can_compute] ** 2
    mfl_flag[can_compute] = FLAG_RHOU2

    # any momentum still missing & short gap -> interpolate
    mfl, m_interp = interp_short_gaps(t, mfl, SHORT_GAP_H)
    mfl_flag[m_interp & (mfl_flag == FLAG_MISSING)] = FLAG_INTERP
    mfl = np.clip(mfl, 0.0, None)

    # =====================================================================
    # 3 & 4. LATENT & SENSIBLE HEAT FLUX (interp <=3h, else Monin-Obukhov)
    # =====================================================================
    # Monin-Obukhov bulk fluxes (Charusombat 2018) from filled T,q,P, using the
    # Eastern Mooring near-surface temperature as the lake surface temp.  The MO
    # fill is applied wherever the mooring Tw (and the met inputs) are available
    # -- i.e. anywhere Hmo/LEmo are finite -- not just inside the winter windows.
    # MO_METHOD chooses the driver: measured/filled u* ("ustar") or the wind
    # speed at the measurement height ("wind", full bulk -- Charusombat 2018).
    n_tw = int(np.isfinite(Tw_arr).sum())
    cs_kw = dict(cool_skin=cool_ok, Rs=Rsf, Rl=Rlf)
    if MO_METHOD == "wind":
        Hmo, LEmo = compute_mo_fluxes_wind(Uf, Tf, qf, Pf, Tw_C=Tw_arr, **cs_kw)
    elif MO_METHOD == "ustar":
        Hmo, LEmo = compute_mo_fluxes(ust, Tf, qf, Pf, Tw_C=Tw_arr, **cs_kw)
    else:
        raise ValueError(f"GLEN_MO_METHOD/MO_METHOD must be 'ustar' or 'wind', "
                         f"got {MO_METHOD!r}")
    mo_arr = {"sensible_heat_flux": Hmo, "latent_heat_flux": LEmo}
    print(f"MO driver: {MO_METHOD!r}; cool_skin={cool_ok}; mooring Tw available "
          f"at {n_tw}/{n}; H finite at {int(np.isfinite(Hmo).sum())}, "
          f"LE finite at {int(np.isfinite(LEmo).sum())}")

    # ---- stand-alone COARE-from-wind diagnostic fluxes -------------------
    # Momentum and turbulent heat fluxes computed SOLELY from the wind speed via
    # the COARE bulk algorithm (u* solved internally; Charusombat 2018), using
    # a mooring's near-surface Tw for the heat fluxes.  These are written
    # as separate variables, independent of the observed/gap-filled fluxes and of
    # MO_METHOD; physical sign retained (not clipped).  Computed once per
    # mooring (Eastern = primary; Southern = sensitivity check) so the two can
    # be compared directly.
    def coare_wind_diagnostics(Tw_arr_diag, reuse_HLE=None):
        _, tau = ustar_from_wind(Uf, Tf, qf, Pf, Tw_C=Tw_arr_diag)
        if reuse_HLE is not None:
            H, LE = reuse_HLE                  # identical call -> reuse
        else:
            H, LE = compute_mo_fluxes_wind(Uf, Tf, qf, Pf, Tw_C=Tw_arr_diag,
                                           **cs_kw)
        return tau, H, LE

    tau_coare, H_coare, LE_coare = coare_wind_diagnostics(
        Tw_arr, reuse_HLE=(Hmo, LEmo) if MO_METHOD == "wind" else None)
    print(f"COARE-from-wind diagnostics (Eastern Mooring Tw): tau finite at "
          f"{int(np.isfinite(tau_coare).sum())}, "
          f"H at {int(np.isfinite(H_coare).sum())}, "
          f"LE at {int(np.isfinite(LE_coare).sum())}")

    tau_coare_south, H_coare_south, LE_coare_south = coare_wind_diagnostics(
        Tw_arr_south)
    print(f"COARE-from-wind diagnostics (Southern Mooring Tw): tau finite at "
          f"{int(np.isfinite(tau_coare_south).sum())}, "
          f"H at {int(np.isfinite(H_coare_south).sum())}, "
          f"LE at {int(np.isfinite(LE_coare_south).sum())}")

    def fill_heat(raw, label):
        y = raw.copy()
        flag = np.where(np.isfinite(y), FLAG_OBS, FLAG_MISSING).astype(np.int8)

        # short-gap interpolation
        y, m_interp = interp_short_gaps(t, y, SHORT_GAP_H)
        flag[m_interp] = FLAG_INTERP

        # Monin-Obukhov bulk fill wherever Tw (hence the MO flux) is available
        mo = mo_arr[label]
        mo_need = (~np.isfinite(y)) & np.isfinite(mo)
        y[mo_need] = mo[mo_need]
        flag[mo_need] = FLAG_MO

        y = np.clip(y, 0.0, None)            # gaps with no Tw remain NaN
        return y, flag

    print("Filling heat fluxes with Monin-Obukhov bulk transfer:")
    le, le_flag = fill_heat(le_raw, "latent_heat_flux")
    sh, sh_flag = fill_heat(sh_raw, "sensible_heat_flux")

    # =====================================================================
    # SUMMARY
    # =====================================================================
    def summarize(name, raw, filled, flag):
        miss0 = int(np.isnan(raw).sum())
        miss1 = int(np.isnan(filled).sum())
        counts = {c: int((flag == c).sum()) for c in
                  (FLAG_OBS, FLAG_INTERP, FLAG_RHOU2, FLAG_WINDBIN, FLAG_MO,
                   FLAG_COARE, FLAG_CLIM, FLAG_MAGNUS, FLAG_PCHIP, FLAG_VECTOR,
                   FLAG_MISSING)}
        print(f"  {name:20s} missing {miss0:6d}->{miss1:6d} | "
              f"obs={counts[0]} interp={counts[1]} rhou2={counts[2]} "
              f"windbin={counts[4]} MO={counts[5]} coareU*={counts[6]} "
              f"clim={counts[7]} magnus={counts[8]} pchip={counts[9]} "
              f"vector={counts[10]} stillNaN={counts[-1]}")

    print("\nFill summary (fluxes):")
    summarize("friction_velocity", ust_raw, ust, ust_flag)
    summarize("momentum_flux", mfl_raw, mfl, mfl_flag)
    summarize("latent_heat_flux", le_raw, le, le_flag)
    summarize("sensible_heat_flux", sh_raw, sh, sh_flag)
    print("Fill summary (met drivers):")
    summarize("air_temperature", T, Tf, T_flag)
    summarize("wind_speed", U, Uf, U_flag)
    summarize("wind_direction", WD, WDf, WD_flag)
    summarize("relative_humidity", RH, RHf, RH_flag)
    summarize("specific_humidity", q, qf, q_flag)
    summarize("air_pressure", P, Pf, P_flag)

    # =====================================================================
    # WRITE OUTPUT (copy everything, overwrite filled vars, add flag vars)
    # =====================================================================
    filled_data = {"friction_velocity": ust, "momentum_flux": mfl,
                   "latent_heat_flux": le, "sensible_heat_flux": sh,
                   "air_temperature": Tf, "wind_speed": Uf,
                   "wind_direction": WDf, "relative_humidity": RHf,
                   "specific_humidity": qf, "air_pressure": Pf}
    flag_data = {"friction_velocity": ust_flag, "momentum_flux": mfl_flag,
                 "latent_heat_flux": le_flag, "sensible_heat_flux": sh_flag,
                 "air_temperature": T_flag, "wind_speed": U_flag,
                 "wind_direction": WD_flag, "relative_humidity": RH_flag,
                 "specific_humidity": q_flag, "air_pressure": P_flag}

    out = nc.Dataset(outfile, "w", format=ds.file_format)
    # dimensions
    for dname, dim in ds.dimensions.items():
        out.createDimension(dname, None if dim.isunlimited() else len(dim))
    # global attrs
    for a in ds.ncattrs():
        out.setncattr(a, getattr(ds, a))
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    notes = (f"{stamp}: gap-filled friction_velocity, momentum_flux, "
             f"latent_heat_flux, sensible_heat_flux, and the meteorological "
             f"drivers air_temperature, wind_speed, wind_direction, "
             f"relative_humidity, specific_humidity, air_pressure (observed "
             f"values unaltered; only gaps filled). Met-driver filling: T/RH/q "
             f"short gaps (<= {SHORT_GAP_H}h) linearly time-interpolated; "
             f"wind_speed & wind_direction gaps <=24h PCHIP-interpolated as a "
             f"wind vector (u,v components, so direction wraps cyclically); "
             f"air_pressure gaps <=8 days PCHIP-interpolated and longer gaps "
             f"filled with monthly climatology; specific_humidity recomputed "
             f"from T,RH,P via Magnus relation where missing. "
             f"momentum_flux=rho*u*^2 (rho=P/(Rd*T*(1+0.608q)), Rd={RD}; dry-air "
             f"fallback q=0 where humidity unavailable); "
             f"u* short gaps interpolated, longer gaps from the COARE bulk "
             f"algorithm solved from wind speed (Charusombat et al. 2018; "
             f"stability from winter-window mooring Tw where known, else "
             f"neutral), overall-median fallback where wind is missing; "
             f"heat fluxes short gaps interpolated, then Monin-Obukhov bulk "
             f"transfer (Charusombat et al. 2018; driver={MO_METHOD} "
             f"[{_MO_DRIVER_DESC}], Eastern-Mooring near-surface Tw"
             f"{', COARE cool-skin correction' if cool_ok else ''}) applied "
             f"wherever the mooring Tw is available (not limited to the winter "
             f"windows). "
             f"Heat-flux gaps with no surface temperature for the MO "
             f"calculation are left as NaN. "
             f"Added stand-alone COARE-from-wind diagnostic variables "
             f"(momentum_flux_coare_wind, sensible_heat_flux_coare_wind, "
             f"latent_heat_flux_coare_wind): momentum and turbulent heat fluxes "
             f"computed solely from the wind speed (Charusombat et al. 2018), "
             f"with Eastern-Mooring Tw for the heat fluxes. Also added "
             f"*_coare_wind_southern variants of the same diagnostic using the "
             f"Southern-Mooring Tw instead, as a sensitivity check on the "
             f"choice of mooring surface temperature. "
             f"See *_fill_method variables for per-point provenance.")
    out.setncattr("gapfill_notes", notes)
    out.setncattr("history", notes + "\n" + getattr(ds, "history", ""))

    # variables
    for vname, var in ds.variables.items():
        fill = var._FillValue if "_FillValue" in var.ncattrs() else None
        newv = out.createVariable(vname, var.dtype, var.dimensions, fill_value=fill)
        for a in var.ncattrs():
            if a == "_FillValue":
                continue
            newv.setncattr(a, getattr(var, a))
        if vname in filled_data:
            d = filled_data[vname]
            newv[:] = d
            # widen valid_min/valid_max to the gap-filled range: gap-fills (esp.
            # MO heat fluxes in strong winter convection) can legitimately exceed
            # the observed QC range, and CF-aware readers (xarray / netCDF4
            # auto-mask) would otherwise mask those values back to NaN
            if np.isfinite(d).any():
                if "valid_min" in var.ncattrs():
                    newv.setncattr("valid_min", float(np.nanmin(d)))
                if "valid_max" in var.ncattrs():
                    newv.setncattr("valid_max", float(np.nanmax(d)))
            prev = newv.qc_notes if "qc_notes" in var.ncattrs() else ""
            newv.setncattr("qc_notes", (prev + " | Gap-filled; see "
                           f"{vname}_fill_method.").strip(" |"))
        else:
            newv[:] = var[:]

    # flag variables
    for vname, fl in flag_data.items():
        fv = f"{vname}_fill_method"
        v = out.createVariable(fv, "i1", ds.variables[vname].dimensions,
                               fill_value=FLAG_FILLVALUE)
        v.setncattr("long_name", f"fill method for {vname}")
        v.setncattr("flag_values", FLAG_VALUES)
        v.setncattr("flag_meanings", FLAG_MEANINGS)
        v.setncattr("comment", FLAG_COMMENT)
        v.setncattr("coordinates", "alt lat lon")
        v[:] = fl

    # stand-alone COARE-from-wind diagnostic flux variables, computed once per
    # mooring (Eastern = primary; Southern = sensitivity check on Tsfc choice)
    def coare_comment_for(mooring_name):
        return (
            "Stand-alone COARE/Monin-Obukhov bulk estimate computed solely from "
            "the wind speed at the measurement height (friction velocity solved "
            "internally; Charusombat et al. 2018), using internally gap-filled "
            "air_temperature, specific_humidity and air_pressure"
            + (" and the COARE cool-skin correction" if cool_ok else "")
            + f". Heat fluxes use the {mooring_name}-Mooring near-surface water "
            "temperature as the surface temperature and are NaN outside its "
            "coverage; momentum uses neutral stability where the mooring Tw is "
            "unavailable. Independent of the observed/gap-filled fluxes; "
            "physical sign retained (not clipped).")

    coare_extra = {
        "momentum_flux_coare_wind": (
            tau_coare, "kg m-1 s-2",
            "Momentum flux from COARE bulk algorithm driven by wind speed "
            "(Eastern-Mooring Tw)", coare_comment_for("Eastern")),
        "sensible_heat_flux_coare_wind": (
            H_coare, "W m-2",
            "Sensible heat flux from COARE bulk algorithm driven by wind speed "
            "(Eastern-Mooring Tw)", coare_comment_for("Eastern")),
        "latent_heat_flux_coare_wind": (
            LE_coare, "W m-2",
            "Latent heat flux from COARE bulk algorithm driven by wind speed "
            "(Eastern-Mooring Tw)", coare_comment_for("Eastern")),
        "momentum_flux_coare_wind_southern": (
            tau_coare_south, "kg m-1 s-2",
            "Momentum flux from COARE bulk algorithm driven by wind speed "
            "(Southern-Mooring Tw sensitivity check)",
            coare_comment_for("Southern")),
        "sensible_heat_flux_coare_wind_southern": (
            H_coare_south, "W m-2",
            "Sensible heat flux from COARE bulk algorithm driven by wind speed "
            "(Southern-Mooring Tw sensitivity check)",
            coare_comment_for("Southern")),
        "latent_heat_flux_coare_wind_southern": (
            LE_coare_south, "W m-2",
            "Latent heat flux from COARE bulk algorithm driven by wind speed "
            "(Southern-Mooring Tw sensitivity check)",
            coare_comment_for("Southern")),
    }
    tdim = ds.variables["time"].dimensions
    for vname, (data, units, long_name, comment) in coare_extra.items():
        v = out.createVariable(vname, "f8", tdim, fill_value=np.nan)
        v.setncattr("long_name", long_name)
        v.setncattr("units", units)
        v.setncattr("comment", comment)
        v.setncattr("coordinates", "alt lat lon")
        v.setncattr("grid_mapping", "crs")
        v.setncattr("cell_methods", "time: mean")
        v[:] = data

    out.close()
    ds.close()
    print(f"\nWrote {outfile}")
    return outfile


if __name__ == "__main__":
    main()
