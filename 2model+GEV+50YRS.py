import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # 只为激活 3D
from matplotlib.collections import PolyCollection
import glob, os, json
from datetime import datetime
from scipy.stats import genextreme as gev

# ========= Basic Utilities =========
def cftime_to_datetime(cftime_obj):
    return datetime(cftime_obj.year, cftime_obj.month, cftime_obj.day,
                    getattr(cftime_obj, "hour", 0),
                    getattr(cftime_obj, "minute", 0),
                    getattr(cftime_obj, "second", 0))

def load_and_extract_pr(model_dir, sub_folder, defined_lat, defined_lon,
                        start_date, end_date):
    """
    读取逐日 pr（kg m-2 s-1 或 mm/day），提取指定格点(最近邻)
    时间窗由 start_date/end_date 控制；返回 (time[], pr_mm_day[])
    """
    ssp_dir  = os.path.join(model_dir, sub_folder)
    data_dir = ssp_dir if os.path.isdir(ssp_dir) else model_dir
    file_paths = sorted(glob.glob(os.path.join(data_dir, "*.nc")))

    all_pr, all_time = [], []
    for fp in file_paths:
        with nc.Dataset(fp, "r") as ds:
            var_pr = ds.variables["pr"]
            pr = var_pr[:]
            time = ds.variables["time"][:]
            lat  = ds.variables["lat"][:]
            lon  = ds.variables["lon"][:]

            units = getattr(var_pr, "units", "").lower()
            mv = getattr(var_pr, "missing_value", None)
            fv = getattr(var_pr, "_FillValue", None)

            # time
            t_units = ds.variables["time"].units
            cal = getattr(ds.variables["time"], "calendar", "standard")
            # Keep cftime datetime objects to support non-Gregorian calendars (e.g., 360_day)
            tdates = np.array(nc.num2date(time, units=t_units, calendar=cal))

            # nearest grid point
            lon_arr = np.array(lon, dtype=float)
            if lon_arr.max() > 180 and defined_lon < 0:
                target_lon = defined_lon + 360.0
            else:
                target_lon = float(defined_lon)
            lat_arr = np.array(lat, dtype=float)
            
            # --- 找到 2 个最近纬度格点、2 个最近经度格点做 2×2 插值 ---
            # 1) 最接近的两个纬度索引
            lat_order = np.argsort(np.abs(lat_arr - float(defined_lat)))
            i0, i1 = np.sort(lat_order[:2])
            y0, y1 = lat_arr[i0], lat_arr[i1]
            if y1 == y0:
                wy0 = 1.0
                wy1 = 0.0
            else:
                wy1 = (defined_lat - y0) / (y1 - y0)
                wy0 = 1.0 - wy1

            # 2) 最接近的两个经度索引
            lon_order = np.argsort(np.abs(lon_arr - target_lon))
            j0, j1 = np.sort(lon_order[:2])
            x0, x1 = lon_arr[j0], lon_arr[j1]
            if x1 == x0:
                wx0 = 1.0
                wx1 = 0.0
            else:
                wx1 = (target_lon - x0) / (x1 - x0)
                wx0 = 1.0 - wx1


            # time window using year range to avoid calendar conversion issues
            years = np.array([t.year for t in tdates], dtype=int)
            mask = (years >= start_date.year) & (years <= end_date.year)
            if not np.any(mask):
                continue

            sel_t = tdates[mask]
            sub_pr = pr[mask, :, :]    # [nt, nlat, nlon]
            
            # 提取四个角点 [nt, 2, 2]
            v00 = sub_pr[:, i0, j0]
            v10 = sub_pr[:, i1, j0]
            v01 = sub_pr[:, i0, j1]
            v11 = sub_pr[:, i1, j1]
            
            # 双线性插值：v = Σ w_lat * w_lon * v_ij
            interp_pr = (v00 * wy0 * wx0 +
                         v10 * wy1 * wx0 +
                         v01 * wy0 * wx1 +
                         v11 * wy1 * wx1).astype(float)

            # lack of testing → NaN
            if mv is not None:
                interp_pr = np.where(interp_pr == mv, np.nan, interp_pr)
            if fv is not None:
                interp_pr = np.where(interp_pr == fv, np.nan, interp_pr)

            # 单位换算：kg m-2 s-1 → ×86400 = mm/day；已是 mm/day 则不变
            if "kg m-2 s-1" in units or "kg m^-2 s^-1" in units or "kg m-2 s-1" in units:
                interp_pr = interp_pr * 86400.0

            all_time.extend(sel_t.tolist())
            all_pr.extend(interp_pr.tolist())

    return np.array(all_time), np.array(all_pr, dtype=float)

# ========= 观测数据读取（最近邻/单格点即可，也可视为插值后的产品） =========
def load_obs_daily_pr(obs_nc_path, var_name, lat_name, lon_name,
                      defined_lat, defined_lon,
                      start_date, end_date):
    
    with nc.Dataset(obs_nc_path, "r") as ds:
        pr  = ds.variables[var_name][:]
        time = ds.variables["time"][:]
        lat  = ds.variables[lat_name][:]
        lon  = ds.variables[lon_name][:]

        units = getattr(ds.variables[var_name], "units", "").lower()
        mv = getattr(ds.variables[var_name], "missing_value", None)
        fv = getattr(ds.variables[var_name], "_FillValue", None)

        t_units = ds.variables["time"].units
        cal     = getattr(ds.variables["time"], "calendar", "standard")
        tdates  = nc.num2date(time, units=t_units, calendar=cal)
        tdates  = np.array([cftime_to_datetime(d) for d in tdates])

        # 最近邻
        lon_arr = np.array(lon, dtype=float)
        lat_arr = np.array(lat, dtype=float)
        if lon_arr.max() > 180 and defined_lon < 0:
            target_lon = defined_lon + 360.0
        else:
            target_lon = float(defined_lon)

        lat_idx = int(np.abs(lat_arr - float(defined_lat)).argmin())
        lon_idx = int(np.abs(lon_arr - target_lon).argmin())

        mask = (tdates >= start_date) & (tdates <= end_date)
        if not np.any(mask):
            return np.array([]), np.array([])

        sel_t  = tdates[mask]
        sel_pr = pr[mask, lat_idx, lon_idx].astype(float)

        if mv is not None:
            sel_pr = np.where(sel_pr == mv, np.nan, sel_pr)
        if fv is not None:
            sel_pr = np.where(sel_pr == fv, np.nan, sel_pr)

        if ("kg m-2 s-1" in units) or ("kg m^-2 s^-1" in units):
            sel_pr = sel_pr * 86400.0

        return sel_t, sel_pr


# ========= QM: 基于历史 GCM vs 观测构建分位映射 =========
def build_qm_lookup(gcm_hist, obs_hist, n_q=100):
    """
    gcm_hist, obs_hist: 历史期 daily 序列 (mm/day)
    返回： (g_q, o_q)，可用 apply_qm(x, g_q, o_q) 进行插值映射。
    """
    gcm_hist = np.asarray(gcm_hist, dtype=float)
    obs_hist = np.asarray(obs_hist, dtype=float)

    mask_g = np.isfinite(gcm_hist)
    mask_o = np.isfinite(obs_hist)
    g = np.sort(gcm_hist[mask_g])
    o = np.sort(obs_hist[mask_o])

    if (len(g) < 10) or (len(o) < 10):
        # 数据太少时，不做 QM，直接返回“恒等映射”
        return np.array([0.0, 1.0]), np.array([0.0, 1.0])

    # 使用统一的概率刻度
    p = np.linspace(0.0, 1.0, n_q)
    g_q = np.quantile(g, p)
    o_q = np.quantile(o, p)
    return g_q, o_q


def apply_qm(series, g_q, o_q):
    """
    对任意 GCM 序列应用已构建的 QM 映射：
    series: 原始 GCM daily(mm/day)
    g_q, o_q: build_qm_lookup 生成的查找表
    """
    series = np.asarray(series, dtype=float)
    # 对 NaN 不做映射
    out = np.full_like(series, np.nan, dtype=float)
    mask = np.isfinite(series)
    if np.all(~mask):
        return out
    out[mask] = np.interp(series[mask], g_q, o_q,
                          left=o_q[0], right=o_q[-1])
    return out


def annual_max(time_arr, pr_mm_day, min_valid_frac=0.9):
    """
    Compute annual maxima from daily precipitation series.

    Parameters
    ----------
    time_arr : array-like
        Can be:
        - cftime / datetime objects (with .year)
        - numpy.datetime64 array
        - int daykey array (YYYYMMDD), e.g., 20141231
        - int year array (YYYY), e.g., 2014
    pr_mm_day : 1D array
        Daily precipitation values (mm/day), same length as time_arr
    min_valid_frac : float
        Minimum fraction of finite values within a year to keep that year

    Returns
    -------
    used_years : 1D int array
    amax : 1D float array
    """
    time_arr = np.asarray(time_arr)
    pr_mm_day = np.asarray(pr_mm_day, dtype=float)

    if time_arr.shape[0] != pr_mm_day.shape[0]:
        raise ValueError(f"time_arr and pr_mm_day length mismatch: {len(time_arr)} vs {len(pr_mm_day)}")

    # ---- 1) Build 'years' array robustly ----
    years = None

    # Case A: int arrays (daykey YYYYMMDD or year YYYY)
    if np.issubdtype(time_arr.dtype, np.integer):
        # Heuristic: if values look like YYYYMMDD (>= 10^7), convert to year by //10000
        # Otherwise treat as year directly.
        vmin = int(np.nanmin(time_arr))
        vmax = int(np.nanmax(time_arr))
        if vmax >= 10_000_000:  # ~1e6 threshold, safe guard
            # likely YYYYMMDD
            years = (time_arr // 10000).astype(int)
        else:
            # likely YYYY already
            years = time_arr.astype(int)

    # Case B: numpy.datetime64
    elif np.issubdtype(time_arr.dtype, np.datetime64):
        years = (time_arr.astype("datetime64[Y]").astype(int) + 1970).astype(int)

    # Case C: object array (cftime / datetime)
    else:
        # If objects have attribute 'year', use it.
        try:
            years = np.array([int(getattr(t, "year")) for t in time_arr], dtype=int)
        except Exception as e:
            raise TypeError("Unsupported time_arr type for annual_max().") from e

    # ---- 2) Compute AMAX per year with completeness check ----
    uy = np.unique(years)
    amax = []
    used_years = []

    for y in uy:
        mask = (years == y)
        vals = pr_mm_day[mask]
        n_all = mask.sum()
        if n_all == 0:
            continue

        valid = np.isfinite(vals)
        if valid.sum() / float(n_all) < float(min_valid_frac):
            continue

        amax.append(np.nanmax(vals))
        used_years.append(int(y))

    amax = np.array(amax, dtype=float)
    used_years = np.array(used_years, dtype=int)

    if len(amax) < 3:
        raise RuntimeError(f"lack of valid years ({len(amax)})")

    return used_years, amax


def fit_gev(amax):
    # SciPy GEV: shape=c=-xi
    c, loc, scale = gev.fit(amax)
    return c, loc, scale

def rl_T_daily(c, loc, scale, T):
    p = 1.0 - 1.0/float(T)
    return float(gev.ppf(p, c, loc=loc, scale=scale))  # mm/day

# ========= Guangzhou Intensity Formulas =========
# q 单位：L/s·ha；t：min；P：年。换算：1 (L/s·ha) = 0.36 (mm/h)
def gz_intensity_Lpsha(P_year, t_min):
    return 13290.63 * (1 + 0.607 * np.log10(P_year)) / ((t_min + 39.126) ** 0.956)

def gz_intensity_mm_per_h(P_year, t_min):
    return gz_intensity_Lpsha(P_year, t_min) * 0.36  # → mm/h

def make_chicago_series(total_depth_mm, duration_min, dt_min=1, r=0.4):
    """
    Generate a Chicago design storm hyetograph.
    total_depth_mm : total storm depth (mm)
    duration_min   : total storm duration (minutes)
    dt_min         : time step (minutes)
    r              : ratio of time to peak (0.4 = 40% before peak)
    """
    
    n = int(duration_min/dt_min)
    tp = int(r * n)
    t = np.arange(0, duration_min+dt_min, dt_min)
    # synthetic symmetrical shape
    a = np.linspace(0, 1, tp, endpoint=False)
    b = np.linspace(1, 0, n - tp + 1)
    pattern = np.concatenate([a, b])
    pattern = pattern / pattern.sum()     # normalize
    intensities = pattern * total_depth_mm / (dt_min/60)  # mm/h
    return t, intensities


def write_timeseries_for_swmm(outfile, name, t, intensities, *, blank_date=True):
    """
    Create a [TIMESERIES] section for SWMM.
    If blank_date=True, leave the Date column blank as requested.
    If fixed_name is provided, use that literal in the Name column (e.g., '10').
    """
    with open(outfile, "w", encoding="utf-8") as f:
        f.write("[TIMESERIES]\n")
        f.write(";;Name           Date       Time       Value     \n")
        f.write(";;-------------- ---------- ---------- ----------\n")
        series_name = str(name)
        for ti, val in zip(t, intensities):
            hh = int(ti // 60)
            mm = int(ti % 60)
            # Date column left blank if requested; preserve Name as given
            if blank_date:
                # spacing mirrors the example: Name, spaces for Date, Time, Value
                f.write(f"{series_name:<14}          {hh}:{mm:02d}       {val:.3f}     \n")
            else:
                f.write(f"{series_name:<14}  01/01/00  {hh}:{mm:02d}  {val:.3f}\n")

'''
def verify_liu6_models_CF(input_dir, output_dir, sub_future="ssp585",
                          hist_range=("1955-01-01","2014-12-31"),
                          fut_range =("2015-01-01","2045-01-01"),
                          lat=23.0, lon=113.0,
                          T=10, duration_hours=2):
    """
    六模型逐一：历史与未来分别 GEV→RL_T，计算 CF_T = RL_T(fut)/RL_T(hist)。
    用广州 IDF 公式算历史 10y–2h，再按 CF_T 缩放到未来；并输出 6 模型中位 CF 的 ensemble 结果。
    """
    os.makedirs(output_dir, exist_ok=True)

    # 6 GCMs according to Liu et al. (2021)
    liu6 = [
        "CanESM5",
        "EC-Earth3",
        "IPSL-CM6A-LR",
        "MPI-ESM1-2-HR",
        "MRI-ESM2-0",
        "NorESM2-MM",
    ]

    # scan available model dirs
    available_dirs = {e.name: e.path for e in os.scandir(input_dir) if e.is_dir()}
    model_paths = {m: available_dirs[m] for m in liu6 if m in available_dirs}
    missing = [m for m in liu6 if m not in model_paths]
    if missing:
        print("⚠️ no model directories found for: ", ", ".join(missing))

    # time ranges
    hstart = datetime.strptime(hist_range[0], "%Y-%m-%d")
    hend   = datetime.strptime(hist_range[1], "%Y-%m-%d")
    fstart = datetime.strptime(fut_range[0], "%Y-%m-%d")
    fend   = datetime.strptime(fut_range[1], "%Y-%m-%d")

    # historical 10y–2h intensity & depth
    tmin = duration_hours * 60.0
    i_hist_mmph = gz_intensity_mm_per_h(T, tmin)
    P_hist_mm   = i_hist_mmph * duration_hours

    summary = {}
    cf_list = []

    for std_name in liu6:
        if std_name not in model_paths:
            continue

        mdir = model_paths[std_name]
        print("\n==============================")
        print(f"Model: {std_name} | Scenario: {sub_future}")
        print("==============================")

        # read daily pr (historical & future)
        t_hist, pr_hist = load_and_extract_pr(mdir, "historical", lat, lon, hstart, hend)
        t_fut , pr_fut  = load_and_extract_pr(mdir, sub_future,  lat, lon, fstart, fend)
        if (len(t_hist)==0) or (len(t_fut)==0):
            print("  !! no data for hist or fut, skipping.")
            continue

        # annual max daily precipitation
        yH, aH = annual_max(t_hist, pr_hist, min_valid_frac=0.9)
        yF, aF = annual_max(t_fut , pr_fut , min_valid_frac=0.9)

        # GEV fitting & RL_T
        cH, locH, scH = fit_gev(aH)
        cF, locF, scF = fit_gev(aF)
        RL_hist = rl_T_daily(cH, locH, scH, T=T)  # mm/day
        RL_fut  = rl_T_daily(cF, locF, scF, T=T)  # mm/day

        # change factor
        CF_T = RL_fut / RL_hist
        i_fut_mmph = i_hist_mmph * CF_T
        P_fut_mm   = P_hist_mm   * CF_T
        cf_list.append(CF_T)

        summary[std_name] = {
            "years_hist": [int(yH.min()), int(yH.max())],
            "years_fut":  [int(yF.min()), int(yF.max())],
            "RL_daily_hist_mm_T": float(RL_hist),
            "RL_daily_fut_mm_T":  float(RL_fut),
            "CF_T": float(CF_T),
            "i_hist_mmph_10y2h": float(i_hist_mmph),
            "i_fut_mmph_10y2h_byCF": float(i_fut_mmph),
            "P_hist_mm_10y2h": float(P_hist_mm),
            "P_fut_mm_10y2h_byCF": float(P_fut_mm),
        }

        print(f"  RL{T} hist/fut = {RL_hist:.2f}/{RL_fut:.2f} mm/day → CF_T={CF_T:.3f}  | "
              f"10y–2h hist {P_hist_mm:.1f} → fut {P_fut_mm:.1f} mm")

    # 6 median CF → ensemble-median future 10y–2h depth
    if cf_list:
        CF_med = float(np.median(cf_list))
        i_fut_med = i_hist_mmph * CF_med
        P_fut_med = P_hist_mm   * CF_med
        summary["_ensemble_median_"] = {
            "CF_T_median_over_6models": CF_med,
            "i_hist_mmph_10y2h": float(i_hist_mmph),
            "i_fut_mmph_10y2h_byCFmedian": float(i_fut_med),
            "P_hist_mm_10y2h": float(P_hist_mm),
            "P_fut_mm_10y2h_byCFmedian": float(P_fut_med),
        }
        print("\n====== Ensemble (6-model) ======")
        print(f"  CF_T median = {CF_med:.3f}  | 10y–2h hist {P_hist_mm:.1f} → fut(median) {P_fut_med:.1f} mm")
    else:
        print("\n⚠️ no valid models processed, skipping ensemble median calculation.")

    # pull the ensemble-median future 10y–2h depth
    P_fut_mm = summary["_ensemble_median_"]["P_fut_mm_10y2h_byCFmedian"]

    # build 2h @ 1-min Chicago series and write TIMESERIES file
    t, i = make_chicago_series(total_depth_mm=P_fut_mm, duration_min=120, dt_min=1, r=0.4)
    out_txt = os.path.join(output_dir, "SSP585_Ensemble10y2h_1min_TIMESERIES.txt")
    write_timeseries_for_swmm(out_txt, "SSP585_Ensemble10y2h", t, i)
    print("SWMM [TIMESERIES] written:", out_txt)

    # save summary JSON
    out_json = os.path.join(output_dir, f"liu6_CF_T_and_10y2h_gz_{sub_future}.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(f"\nSaved summary to: {out_json}")
'''
    
    
def verify_liu6_models_CF_ensemble_qm(
        input_dir, output_dir,
        #obs_nc_path,           # 观测 NetCDF 路径（你需要自己改）
        #obs_var="pr",          # 观测变量名
        #obs_lat_name="lat",    # 观测文件的纬度变量名
        #obs_lon_name="lon",    # 观测文件的经度变量名
        sub_future="ssp585",
        hist_range=("1955-01-01","2014-12-31"),
        fut_range =("2015-01-01","2045-01-01"),
        lat=23.0, lon=113.0,
        T=10, duration_hours=2):
    """
    实现你要求的：
    (1) 读取6个 GCM daily rainfall (historical + future)，做 2×2 双线性插值；
    (2) 用观测做 QM 偏差订正；
    (3) 在 daily 序列上做 ensemble median（历史、未来各一条）；
    (4) 对 ensemble 序列做 AMAX→GEV→CF；
    (5) 输出 ensemble 的 Chicago 10y–2h 时间序列，
        并额外输出每个 GCM 自己的 Chicago 10y–2h 时间序列（已 QM）。
    """
    def time_to_daykey_arr(time_arr):
        """
        Convert cftime/datetime-like array to int day keys: YYYYMMDD.
        This avoids comparing mixed cftime calendars in numpy set ops.
        """
        keys = np.empty(len(time_arr), dtype=np.int32)
        for i, t in enumerate(time_arr):
            y = int(getattr(t, "year"))
            m = int(getattr(t, "month"))
            d = int(getattr(t, "day"))
            keys[i] = y * 10000 + m * 100 + d
        return keys

    def align_series_by_daykey(pr_arr, daykey_arr, common_daykeys):
        """
        Align a 1D precipitation series to common_daykeys order.
        Ensures same length and same day order across all GCMs.
        """
        # Build key -> index mapping (if duplicates exist, keep the first occurrence)
        idx_map = {}
        for j, k in enumerate(daykey_arr):
            kk = int(k)
            if kk not in idx_map:
                idx_map[kk] = j

        aligned = np.empty(len(common_daykeys), dtype=float)
        for j, k in enumerate(common_daykeys):
            kk = int(k)
            aligned[j] = pr_arr[idx_map[kk]]
        return aligned

    
    future_amax_dict = {}  
    os.makedirs(output_dir, exist_ok=True)

    # 动态扫描 input_dir 下的所有模型子文件夹
    available_dirs = {e.name: e.path for e in os.scandir(input_dir) if e.is_dir()}
    # 排除明显的非模型文件夹（可按需扩展）
    exclude = {"historical", "ssp585", "ssp245", "docs", "output"}
    all_models = sorted([name for name in available_dirs.keys() if name not in exclude])
    model_paths = {m: available_dirs[m] for m in all_models}
    if not all_models:
        print("⚠️ 未在 input_dir 中发现任何模型子文件夹。")

    # 时间
    hstart = datetime.strptime(hist_range[0], "%Y-%m-%d")
    hend   = datetime.strptime(hist_range[1], "%Y-%m-%d")
    fstart = datetime.strptime(fut_range[0], "%Y-%m-%d")
    fend   = datetime.strptime(fut_range[1], "%Y-%m-%d")

    # 历史广州 10y–2h 基线
    tmin = duration_hours * 60.0
    i_hist_mmph = gz_intensity_mm_per_h(T, tmin)
    P_hist_mm   = i_hist_mmph * duration_hours

    '''
    # 读取观测历史
    print("\n=== Load OBS daily series for QM ===")
    t_obs_hist, pr_obs_hist = load_obs_daily_pr(
        obs_nc_path, obs_var, obs_lat_name, obs_lon_name,
        lat, lon, hstart, hend
    )
    if len(t_obs_hist) == 0:
        print("⚠️ 观测历史期为空，QM 将退化为恒等映射。")
    '''

    # 存放各 GCM 经 QM 订正后的 daily 序列
    hist_series_bc = []  # list of 1D arrays
    fut_series_bc  = []
    hist_times     = []  # 新增：存每个 GCM 的历史 time 数组
    fut_times      = []  # 新增：存每个 GCM 的未来 time 数组
    #t_hist_ref = None
    #t_fut_ref  = None

    summary = {}
    cf_list = []

    for std_name in all_models:
        if std_name not in model_paths:
            continue
        mdir = model_paths[std_name]
        print("\n==============================")
        print(f"Model: {std_name} | Scenario: {sub_future}")
        print("==============================")

        # 1) 读取 GCM 历史/未来 daily（双线性插值后）
        t_hist_raw, pr_hist_raw = load_and_extract_pr(
            mdir, "historical", lat, lon, hstart, hend
        )
        t_fut_raw, pr_fut_raw = load_and_extract_pr(
            mdir, sub_future, lat, lon, fstart, fend
        )
        if (len(t_hist_raw) == 0) or (len(t_fut_raw) == 0):
            print("  !! 历史或未来数据缺失，跳过该 GCM")
            continue
        
                # 记录 time 轴，后面统一求交集做对齐
        hist_times.append(time_to_daykey_arr(t_hist_raw))
        fut_times.append(time_to_daykey_arr(t_fut_raw))


        '''
        # 对齐时间轴（NEX-GDDP 一般是整齐 daily），这里简单假定是完整日序列
        if t_hist_ref is None:
            t_hist_ref = t_hist_raw
        else:
            if not np.array_equal(t_hist_ref, t_hist_raw):
                print("  ⚠️ 历史 time 轴不一致，简单截取重叠部分。")
        if t_fut_ref is None:
            t_fut_ref = t_fut_raw
        else:
            if not np.array_equal(t_fut_ref, t_fut_raw):
                print("  ⚠️ 未来 time 轴不一致，简单截取重叠部分。")
        
        # 2) 构建 QM 映射（GCM_hist vs OBS_hist），然后作用到 GCM_hist & GCM_fut
        if len(pr_obs_hist) > 0:
            g_q, o_q = build_qm_lookup(pr_hist_raw, pr_obs_hist)
            pr_hist_bc = apply_qm(pr_hist_raw, g_q, o_q)
            pr_fut_bc  = apply_qm(pr_fut_raw,  g_q, o_q)
        else:
            # 没有观测 → 不做 QM
            pr_hist_bc = pr_hist_raw.copy()
            pr_fut_bc  = pr_fut_raw.copy()

        hist_series_bc.append(pr_hist_bc)
        fut_series_bc.append(pr_fut_bc)
        '''
        # 存入 ensemble 容器（未做 QM）
        hist_series_bc.append(pr_hist_raw)
        fut_series_bc.append(pr_fut_raw)

        # 3) 对单个 GCM（已 QM）做 AMAX → GEV → RL → CF，并输出单独的 Chicago 序列
        yH, aH = annual_max(t_hist_raw, pr_hist_raw, min_valid_frac=0.9)
        yF, aF = annual_max(t_fut_raw,  pr_fut_raw,  min_valid_frac=0.9)
        future_amax_dict[std_name] = (yF, aF)


        cH, locH, scH = fit_gev(aH)
        cF, locF, scF = fit_gev(aF)
        RL_hist = rl_T_daily(cH, locH, scH, T=T)
        RL_fut  = rl_T_daily(cF, locF, scF, T=T)
        CF_T    = RL_fut / RL_hist
        cf_list.append(CF_T)

        i_fut_mmph = i_hist_mmph * CF_T
        P_fut_mm   = P_hist_mm   * CF_T

        summary[std_name] = {
            "years_hist": [int(yH.min()), int(yH.max())],
            "years_fut":  [int(yF.min()), int(yF.max())],
            "RL_daily_hist_mm_T": float(RL_hist),
            "RL_daily_fut_mm_T":  float(RL_fut),
            "CF_T": float(CF_T),
            "i_hist_mmph_10y2h": float(i_hist_mmph),
            "i_fut_mmph_10y2h_byCF": float(i_fut_mmph),
            "P_hist_mm_10y2h": float(P_hist_mm),
            "P_fut_mm_10y2h_byCF": float(P_fut_mm),
        }

        print(f"  [Single GCM] RL{T} hist/fut = {RL_hist:.2f}/{RL_fut:.2f} mm/day "
              f"→ CF_T={CF_T:.3f} | 10y–2h hist {P_hist_mm:.1f} → fut {P_fut_mm:.1f} mm")

        # 生成该 GCM 自己的 Chicago 10y–2h 时间序列
        t_chi, i_chi = make_chicago_series(
            total_depth_mm=P_fut_mm, duration_min=int(duration_hours*60),
            dt_min=1, r=0.4
        )
        out_txt_model = os.path.join(
            output_dir, f"{std_name}_{sub_future}_10y2h_1min_TIMESERIES.txt"
        )
        write_timeseries_for_swmm(
            out_txt_model,
            f"{std_name}_10y2h",
            t_chi, i_chi,
            blank_date=True
        )
        print("  SWMM [TIMESERIES] written for single GCM:", out_txt_model)

    # ========== 在 daily 序列层面做 ensemble median ==========
    if (not hist_series_bc) or (not fut_series_bc):
        print("⚠️ 无可用 GCM 序列，无法生成 ensemble。")
        return

    # ========== 在 daily 序列层面做 ensemble median（按 daykey 对齐） ==========
    if (not hist_series_bc) or (not fut_series_bc):
        print("⚠️ 无可用 GCM 序列，无法生成 ensemble。")
        return

    # 1) 求“所有 GCM 共同拥有的历史日期(daykey)”交集
    common_hist_key = hist_times[0]
    for th_key in hist_times[1:]:
        common_hist_key = np.intersect1d(common_hist_key, th_key)

    # 2) 求“所有 GCM 共同拥有的未来日期(daykey)”交集
    common_fut_key = fut_times[0]
    for tf_key in fut_times[1:]:
        common_fut_key = np.intersect1d(common_fut_key, tf_key)

    # 重要：排序，保证所有模型按同一日期顺序对齐（intersect1d 一般会返回排序结果，但这里明确一下）
    common_hist_key = np.sort(common_hist_key)
    common_fut_key  = np.sort(common_fut_key)

    # 3) 把每个 GCM 的序列对齐到共同日期（严格按 common_*_key 的顺序）
    aligned_hist = []
    aligned_fut  = []
    for th_key, sh, tf_key, sf in zip(hist_times, hist_series_bc, fut_times, fut_series_bc):
        aligned_hist.append(align_series_by_daykey(sh, th_key, common_hist_key))
        aligned_fut.append(align_series_by_daykey(sf, tf_key, common_fut_key))

    # 4) 现在所有数组长度一致，可以安全 vstack
    hist_stack = np.vstack(aligned_hist)  # [n_models, n_days_hist_common]
    fut_stack  = np.vstack(aligned_fut)   # [n_models, n_days_fut_common]

    # 交集 daykey 轴作为 ensemble 的时间轴（用于 annual_max 按年份分组）
    t_hist_ref = common_hist_key
    t_fut_ref  = common_fut_key


    pr_hist_ens = np.median(hist_stack, axis=0)
    pr_fut_ens  = np.median(fut_stack,  axis=0)

    # 5) 对 ensemble daily 序列做 AMAX → GEV → CF
    yH_e, aH_e = annual_max(t_hist_ref, pr_hist_ens, min_valid_frac=0.9)
    yF_e, aF_e = annual_max(t_fut_ref,  pr_fut_ens,  min_valid_frac=0.9)

    cH_e, locH_e, scH_e = fit_gev(aH_e)
    cF_e, locF_e, scF_e = fit_gev(aF_e)

    RL_hist_ens = rl_T_daily(cH_e, locH_e, scH_e, T=T)
    RL_fut_ens  = rl_T_daily(cF_e, locF_e, scF_e, T=T)
    CF_ens      = RL_fut_ens / RL_hist_ens

    i_fut_mmph_ens = i_hist_mmph * CF_ens
    P_fut_mm_ens   = P_hist_mm   * CF_ens

    summary["_ensemble_field_median_"] = {
        "RL_daily_hist_mm_T": float(RL_hist_ens),
        "RL_daily_fut_mm_T":  float(RL_fut_ens),
        "CF_T_ensemble_field": float(CF_ens),
        "i_hist_mmph_10y2h": float(i_hist_mmph),
        "i_fut_mmph_10y2h_byCF": float(i_fut_mmph_ens),
        "P_hist_mm_10y2h": float(P_hist_mm),
        "P_fut_mm_10y2h_byCF": float(P_fut_mm_ens),
    }

    print("\n====== Ensemble (field-level median over 6 GCMs) ======")
    print(f"  RL{T} hist/fut = {RL_hist_ens:.2f}/{RL_fut_ens:.2f} mm/day "
          f"→ CF_ens={CF_ens:.3f} | 10y–2h hist {P_hist_mm:.1f} → fut(ens) {P_fut_mm_ens:.1f} mm")

    # 6) ensemble Chicago 10y–2h 时间序列
    t_ens, i_ens = make_chicago_series(
        total_depth_mm=P_fut_mm_ens, duration_min=int(duration_hours*60),
        dt_min=1, r=0.4
    )
    out_txt_ens = os.path.join(output_dir, f"{sub_future}_Ensemble10y2h_1min_TIMESERIES_QM.txt")
    write_timeseries_for_swmm(
        out_txt_ens,
        f"{sub_future}_Ens10y2h",
        t_ens, i_ens,
        blank_date=True
    )
    print("SWMM [TIMESERIES] written for ensemble:", out_txt_ens)

    # 保存汇总 JSON
    out_json = os.path.join(output_dir, f"liu6_CF_T_and_10y2h_gz_{sub_future}_ensemble_QM.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print("\nSaved summary to:", out_json)
    out_png = os.path.join(output_dir, f"Future_AMAX_3D_ssp585.png")
    plot_future_amax_3d(future_amax_dict, out_png)


def plot_hist_vs_future_amax(input_dir, model_name,
                             lat=23.0, lon=113.0,
                             hist_range=("1955-01-01","2014-12-31"),
                             fut_range =("2015-01-01","2045-01-01"),
                             show_year_window=(1955, 2045),
                             out_dir=".", tag="SSP585"):
    """
    仅区分历史与未来：绘制年最大日雨（AMAX）的柱状/折线图，
    并分别给出历史与未来的均值水平线。
    需要你已有的 load_and_extract_pr() 与 annual_max()。
    """
    os.makedirs(out_dir, exist_ok=True)
    hstart = datetime.strptime(hist_range[0], "%Y-%m-%d")
    hend   = datetime.strptime(hist_range[1], "%Y-%m-%d")
    fstart = datetime.strptime(fut_range[0], "%Y-%m-%d")
    fend   = datetime.strptime(fut_range[1], "%Y-%m-%d")

    mdir = os.path.join(input_dir, model_name)

    # read daily pr (historical & future)
    t_hist, pr_hist = load_and_extract_pr(mdir, "historical", lat, lon, hstart, hend)
    t_fut , pr_fut  = load_and_extract_pr(mdir, "ssp585",   lat, lon, fstart, fend)

    # annual max daily precipitation
    yH, aH = annual_max(t_hist, pr_hist, min_valid_frac=0.9)
    yF, aF = annual_max(t_fut , pr_fut ,  min_valid_frac=0.9)

    
    y0, y1 = show_year_window
    maskH = (yH >= y0) & (yH <= y1)
    maskF = (yF >= y0) & (yF <= y1)
    yH_s, aH_s = yH[maskH], aH[maskH]
    yF_s, aF_s = yF[maskF], aF[maskF]

    plt.figure(figsize=(12, 4.2))
    # historical
    plt.bar(yH_s, aH_s, width=0.8, alpha=0.35, edgecolor='none', label='Historical AMAX')
    plt.plot(yH_s, aH_s, lw=1.2)
    # future
    plt.bar(yF_s, aF_s, width=0.8, alpha=0.35, edgecolor='none', label=f'Future {tag} AMAX')
    plt.plot(yF_s, aF_s, lw=1.2)

    # average lines
    if len(aH_s) > 0:
        meanH = float(np.mean(aH_s))
        plt.hlines(meanH, yH_s.min()-0.4, yH_s.max()+0.4, linestyles='dashed', linewidth=1.8)
        plt.text((yH_s.min()+yH_s.max())/2, meanH*1.02, f"Hist mean: {meanH:.1f} mm",
                 ha='center', va='bottom', fontsize=9)
    if len(aF_s) > 0:
        meanF = float(np.mean(aF_s))
        plt.hlines(meanF, yF_s.min()-0.4, yF_s.max()+0.4, linestyles='dashed', linewidth=1.8)
        plt.text((yF_s.min()+yF_s.max())/2, meanF*1.02, f"Future mean: {meanF:.1f} mm",
                 ha='center', va='bottom', fontsize=9)

    # historical/future split line (using future start year)
    split_year = datetime.strptime(fut_range[0], "%Y-%m-%d").year
    plt.axvline(split_year, color='gray', linestyle='--', linewidth=1.2)
    plt.text(split_year+0.2, plt.ylim()[1]*0.95, f"{split_year}", rotation=90, va='top', fontsize=9)

    plt.xlim(y0-1, y1+1)
    plt.xlabel("Year")
    plt.ylabel("Annual maximum daily precipitation (mm)")
    plt.title(f"{model_name} — AMAX (Historical vs {tag})")
    plt.grid(alpha=0.25)
    plt.legend(loc='upper left', frameon=False)
    plt.tight_layout()

    out_png = os.path.join(out_dir, f"{model_name}_AMAX_hist_vs_{tag}.png")
    plt.savefig(out_png, dpi=200)
    plt.close()
    print("Saved:", out_png)
    
def plot_future_amax_3d(future_amax_dict, out_png_path):
    """
    画出 6 个 GCM 的未来 AMAX 曲线的 3D “层叠图”，并保存为 PNG。

    future_amax_dict: {model_name: (years_array, amax_array)}
    out_png_path: PNG 文件完整路径
    """
    # 保证模型顺序固定
    model_names = list(future_amax_dict.keys())

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection="3d")

    polys = []
    ys = []   # 3D 里的“Y 轴位置”（模型索引）

    for j, name in enumerate(model_names):
        years, amax = future_amax_dict[name]

        # 按年份排序，避免折线乱跳
        order = np.argsort(years)
        xs = years[order]
        zs = amax[order]

        # 构造 2D 顶点 (x,z)，再封闭到底边，形成一个“填充区域”
        verts = list(zip(xs, zs))
        verts.append((xs[-1], 0.0))
        verts.append((xs[0], 0.0))

        polys.append(verts)
        ys.append(j)  # 每个模型放在 y=j 的平面上

    poly = PolyCollection(polys, alpha=0.6)
    # 这里不特别指定颜色，用默认色系即可；如果想统一橙色可以加 facecolor='orange'
    ax.add_collection3d(poly, zs=ys, zdir='y')

    ax.set_xlabel("Year")
    ax.set_ylabel("GCM index")
    ax.set_zlabel("Annual maximum daily precipitation (mm)")

    ax.set_yticks(range(len(model_names)))
    ax.set_yticklabels(model_names)

    # x 轴范围统一一下（根据所有年份）
    all_years = np.concatenate([future_amax_dict[name][0] for name in model_names])
    ax.set_xlim(all_years.min(), all_years.max())

    # z 轴从 0 开始更好看
    all_amax = np.concatenate([future_amax_dict[name][1] for name in model_names])
    ax.set_zlim(0, all_amax.max() * 1.1)

    ax.view_init(elev=25, azim=-60)  # 调一下视角

    plt.tight_layout()
    fig.savefig(out_png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("3D future AMAX figure saved to:", out_png_path)

def main():
    
    input_dir  = r"C:\Users\yxz2315\Desktop\pythonvscode\model_downloads"
    output_dir = r"C:\Users\yxz2315\Desktop\pythonvscode\model_outputs"

    # coordinates for Guangzhou
    lat, lon = 23.0, 113.0

    # historical 60 years and future 30 years (consistent with the paper)
    hist_range = ("1955-01-01", "2014-12-31")
    fut_range  = ("2015-01-01", "2045-01-01")

    # target event: 10-year return period, 2 hours
    T = 10
    duration_hours = 2

    verify_liu6_models_CF_ensemble_qm(
        input_dir=input_dir, output_dir=output_dir,
        sub_future="ssp585",
        hist_range=hist_range, fut_range=fut_range,
        lat=lat, lon=lon,
        T=T, duration_hours=duration_hours
    )
    plot_hist_vs_future_amax(
        input_dir, "NorESM2-MM",
        lat=23.0, lon=113.0,
        hist_range=("1955-01-01","2014-12-31"),
        fut_range =("2015-01-01","2045-01-01"),
        show_year_window=(1955, 2045),
        out_dir=output_dir, tag="SSP585"
    )
    # 1) 指定输出目录和 JSON 文件名
    output_dir = r"C:\Users\yxz2315\Desktop\pythonvscode\model_outputs"
    json_name  = "liu6_CF_T_and_10y2h_gz_ssp585_ensemble_QM.json"  # 或 _noQM.json
    json_path  = os.path.join(output_dir, json_name)

    # 2) 读取 summary
    with open(json_path, "r", encoding="utf-8") as f:
        summary = json.load(f)

    # 3) 从 summary 动态获取模型名称（排除 ensemble 键）
    ensemble_keys = {
        "_ensemble_field_median_",
        "_ensemble_field_median_noQM_",
        "_ensemble_median_",
    }
    models = sorted([k for k in summary.keys() if k not in ensemble_keys])

    # 4) 取一个历史基准（对所有模型一样，取第一个有的）
    P_hist = None
    for m in models:
        if m in summary:
            P_hist = summary[m]["P_hist_mm_10y2h"]
            break

    if P_hist is None:
        raise RuntimeError("找不到任何模型的 P_hist_mm_10y2h，检查 JSON。")

    # 5) 取 ensemble future 10y–2h
    #   对 QM 版本：key 通常是 "_ensemble_field_median_"
    #   对 no-QM 版本：key 通常是 "_ensemble_field_median_noQM_"
    if "_ensemble_field_median_" in summary:
        ens_key = "_ensemble_field_median_"
    elif "_ensemble_field_median_noQM_" in summary:
        ens_key = "_ensemble_field_median_noQM_"
    elif "_ensemble_median_" in summary:
        ens_key = "_ensemble_median_"
    else:
        raise RuntimeError("找不到 ensemble 汇总 key，检查 summary 结构。")

    P_fut_ens = summary[ens_key]["P_fut_mm_10y2h_byCF"]

    # 6) 取各模型各自的未来 10y–2h
    P_fut_list = []
    available_models = []
    for m in models:
        if m in summary:
            P_fut_list.append(summary[m]["P_fut_mm_10y2h_byCF"])
            available_models.append(m)
        else:
            print(f"⚠️ 警告: {m} 不在 JSON 里，跳过。")

    # 7) 组合成 8 根柱子的标签和数值
    labels = ["Hist-10y2h", "Future-Ens"] + available_models
    values = [P_hist, P_fut_ens] + P_fut_list

    # 8) 画柱状图
    x = np.arange(len(labels))

    plt.figure(figsize=(10, 5))
    plt.bar(x, values)
    plt.xticks(x, labels, rotation=30, ha="right")
    plt.ylabel("Total 2h rainfall (mm)")
    plt.title("Historical vs future 10-year 2h rainfall (Guangzhou, ssp585)")
    plt.tight_layout()
    plt.grid(axis="y", linestyle="--", alpha=0.3)

    plt.show()
    # save figure as PNG
    out_png = os.path.join(output_dir, "Historical_vs_Future_10y2h_Rainfall_ssp585.png")
    plt.savefig(out_png, dpi=200)
    

if __name__ == "__main__":
    main()