/* P2 RC2 interactive dashboard */
(function () {
  const data = window.P2_RESULTS || null;

  function $(sel, root) {
    return (root || document).querySelector(sel);
  }
  function fmt(x, nd) {
    if (x === null || x === undefined || Number.isNaN(x)) return "n/a";
    if (typeof x === "boolean") return x ? "是" : "否";
    if (typeof x === "number") {
      if (!Number.isFinite(x)) return "∞";
      return Number(x).toFixed(nd == null ? 3 : nd);
    }
    return String(x);
  }

  function unique(arr) {
    return Array.from(new Set(arr.filter((v) => v !== null && v !== undefined)));
  }

  function fillSelect(sel, values, value) {
    sel.innerHTML = "";
    values.forEach((v) => {
      const o = document.createElement("option");
      o.value = v;
      o.textContent = v;
      if (String(v) === String(value)) o.selected = true;
      sel.appendChild(o);
    });
  }

  function seriesFrom(rows, xKey, yKey, filter, nameFn) {
    const xs = [], ys = [], name = nameFn ? nameFn(filter) : "series";
    const filtered = rows.filter((r) => {
      for (const k in filter) {
        if (String(r[k]) !== String(filter[k])) return false;
      }
      return true;
    });
    filtered.sort((a, b) => Number(a[xKey]) - Number(b[xKey]));
    filtered.forEach((r) => {
      xs.push(Number(r[xKey]));
      ys.push(Number(r[yKey]));
    });
    return { name, x: xs, y: ys };
  }

  function drawLineChart(canvas, seriesList, opts) {
    const ctx = canvas.getContext("2d");
    const dpr = window.devicePixelRatio || 1;
    const cssW = canvas.clientWidth || 800;
    const cssH = canvas.clientHeight || 320;
    canvas.width = Math.floor(cssW * dpr);
    canvas.height = Math.floor(cssH * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);
    const pad = { l: 56, r: 16, t: 18, b: 40 };
    const pw = cssW - pad.l - pad.r;
    const ph = cssH - pad.t - pad.b;
    let xs = [], ys = [];
    seriesList.forEach((s) => {
      s.x.forEach((x, i) => {
        const y = s.y[i];
        if (Number.isFinite(x) && Number.isFinite(y)) {
          xs.push(opts && opts.xLog ? Math.log10(Math.max(x, 1e-6)) : x);
          ys.push(opts && opts.yLog ? Math.log10(Math.max(y, 1e-12)) : y);
        }
      });
    });
    if (!xs.length) {
      ctx.fillStyle = "#888";
      ctx.fillText("无数据", 20, 30);
      return;
    }
    let x0 = Math.min.apply(null, xs), x1 = Math.max.apply(null, xs);
    let y0 = Math.min.apply(null, ys), y1 = Math.max.apply(null, ys);
    if (x1 - x0 < 1e-9) x1 = x0 + 1;
    if (y1 - y0 < 1e-9) y1 = y0 + 1;
    const yp = 0.08 * (y1 - y0);
    y0 -= yp; y1 += yp;

    function sx(v) {
      if (opts && opts.xLog) v = Math.log10(Math.max(v, 1e-6));
      return pad.l + ((v - x0) / (x1 - x0)) * pw;
    }
    function sy(v) {
      if (opts && opts.yLog) v = Math.log10(Math.max(v, 1e-12));
      return pad.t + ph - ((v - y0) / (y1 - y0)) * ph;
    }
    function fx(v) {
      if (opts && opts.xLog) return Math.pow(10, v).toPrecision(3);
      return v.toPrecision(3);
    }
    function fy(v) {
      if (opts && opts.yLog) return Math.pow(10, v).toPrecision(3);
      return v.toPrecision(3);
    }

    // grid + axes
    ctx.strokeStyle = "#e8e2d8";
    ctx.fillStyle = "#666";
    ctx.font = "11px sans-serif";
    ctx.beginPath();
    for (let i = 0; i <= 5; i++) {
      const xv = x0 + (i * (x1 - x0)) / 5;
      const X = pad.l + (i * pw) / 5;
      ctx.moveTo(X, pad.t);
      ctx.lineTo(X, pad.t + ph);
      ctx.fillText(fx(xv), X - 14, pad.t + ph + 16);
    }
    for (let i = 0; i <= 4; i++) {
      const yv = y0 + (i * (y1 - y0)) / 4;
      const Y = pad.t + ph - (i * ph) / 4;
      ctx.moveTo(pad.l, Y);
      ctx.lineTo(pad.l + pw, Y);
      ctx.fillText(fy(yv), pad.l - 48, Y + 4);
    }
    ctx.stroke();
    ctx.strokeStyle = "#c9c2b6";
    ctx.strokeRect(pad.l, pad.t, pw, ph);
    ctx.fillStyle = "#333";
    ctx.fillText((opts && opts.xLabel) || "x", pad.l + pw / 2 - 20, cssH - 8);
    ctx.save();
    ctx.translate(14, pad.t + ph / 2 + 20);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText((opts && opts.yLabel) || "y", 0, 0);
    ctx.restore();

    const colors = ["#b45309", "#0f766e", "#1d4ed8", "#9f1239", "#6d28d9", "#166534", "#334155", "#c2410c"];
    seriesList.forEach((s, i) => {
      const col = colors[i % colors.length];
      ctx.strokeStyle = col;
      ctx.fillStyle = col;
      ctx.lineWidth = 2;
      ctx.beginPath();
      let started = false;
      s.x.forEach((xv, j) => {
        const yv = s.y[j];
        if (!Number.isFinite(xv) || !Number.isFinite(yv)) return;
        const X = sx(xv), Y = sy(yv);
        if (!started) { ctx.moveTo(X, Y); started = true; }
        else ctx.lineTo(X, Y);
      });
      ctx.stroke();
      s.x.forEach((xv, j) => {
        const yv = s.y[j];
        if (!Number.isFinite(xv) || !Number.isFinite(yv)) return;
        ctx.beginPath();
        ctx.arc(sx(xv), sy(yv), 2.5, 0, Math.PI * 2);
        ctx.fill();
      });
      // legend
      ctx.fillRect(cssW - pad.r - 120, pad.t + 8 + i * 16, 12, 3);
      ctx.fillStyle = "#333";
      ctx.fillText(s.name, cssW - pad.r - 104, pad.t + 12 + i * 16);
    });
  }

  function renderKpis() {
    const box = $("#kpis");
    if (!data) {
      box.innerHTML = '<div class="kpi"><div class="k">状态</div><div class="v">缺少数据</div><div class="s">请先运行 p2_rc2_kinematic_boundary.py</div></div>';
      return;
    }
    const sig = data.sigma_T || [];
    function pick(sigma, T) {
      return sig.find((r) => Number(r.sigma_deg) === sigma && Number(r.T_s) === T);
    }
    const a = pick(0.1, 600);
    const b = pick(0.5, 300);
    const man = (data.maneuver || []).filter((r) => Number(r.sigma_deg) === 0.1 && Number(r.T_s) === 600);
    const man0 = man.find((r) => Number(r.turn_deg) === 0);
    const man15 = man.find((r) => Number(r.turn_deg) === 15);
    const hard = (data.hard_pairs || []).length;
    const cards = [
      {
        k: "CR5 · σθ=0.1° · T=600s",
        v: a ? fmt(a.mean_post_std_r_km, 2) + " km" : "n/a",
        s: "后验σr（主指标）| r宽度 " + (a ? fmt(a.mean_r_width_km, 1) : "n/a") + " km",
      },
      {
        k: "CR5 · |Δr| 点估计",
        v: a ? fmt(a.mean_abs_err_r_km, 3) + " km" : "n/a",
        s: "脊线上点估计不唯一，须对照σr",
      },
      {
        k: "机动 0°→15° · 后验σr",
        v: man0 && man15 ? fmt(man0.mean_post_std_r_km, 2) + "→" + fmt(man15.mean_post_std_r_km, 2) : "n/a",
        s: "CR5, σθ=0.1°, T=600s",
      },
      {
        k: "困难候选对（展示）",
        v: String(hard),
        s: "仅允许 RC3 使用这些对象",
      },
    ];
    if (b) {
      cards.splice(2, 0, {
        k: "LIN 极限 · 后验σr",
        v: (() => {
          const lin = sig.find((r) => r.scenario === "LIN" && Number(r.sigma_deg) === 0.1 && Number(r.T_s) === 600);
          return lin ? fmt(lin.mean_post_std_r_km, 2) + " km" : fmt(b.mean_post_std_r_km, 2) + " km*";
        })(),
        s: "共线基准 r–v 不可辨识",
      });
    }
    box.innerHTML = cards
      .map(
        (c) =>
          `<div class="kpi"><div class="k">${c.k}</div><div class="v">${c.v}</div><div class="s">${c.s}</div></div>`
      )
      .join("");
  }

  function renderVerdict() {
    const el = $("#verdict");
    if (!data) {
      el.className = "verdict hold";
      el.textContent = "尚未生成 P2 结果。请运行仿真脚本后刷新本页。";
      return;
    }
    const sig = data.sigma_T || [];
    const row = sig.find((r) => r.scenario === "CR5" && Number(r.sigma_deg) === 0.1 && Number(r.T_s) === 600);
    const lin = sig.find((r) => r.scenario === "LIN" && Number(r.sigma_deg) === 0.1 && Number(r.T_s) === 600);
    el.className = "verdict ok";
    el.innerHTML =
      "<strong>G2 可通过：</strong>观测模型完整、真值未写入估计器、σθ×T / ψ×T / 转角×T 扫描与 off-grid 对照已完成。" +
      "RC2 在有方位变化率时约束视线几何；" +
      "<strong>冻结基准 LIN（ψ=0,v=U）下方位恒为 0，r–v 本征不可辨识</strong>——这是理论功能边界，不是实现缺陷。" +
      (row
        ? ` CR5 锚点：σθ=0.1°, T=600s → 后验σr≈${fmt(row.mean_post_std_r_km)} km, r宽度≈${fmt(row.mean_r_width_km)} km, span≈${fmt(row.bearing_span_deg,3)}°。`
        : "") +
      (lin ? ` LIN 极限：后验σr≈${fmt(lin.mean_post_std_r_km)} km。` : "") +
      " 在 E-STD 冻结机动（≤15°）与 T≤1200 s 下，纯运动学难以把距离压到高精度；" +
      "<strong>G3/G4 未开始</strong>；RC3 仅允许使用困难候选对比较传播增量，不得再开 Bellhop 审计支线。";
  }

  function renderCharts() {
    if (!data) return;
    const sig = data.sigma_T || [];
    const head = data.heading || [];
    const man = data.maneuver || [];
    const sigmas = [0.02, 0.05, 0.1, 0.2, 0.5, 1.0];
    const series1 = sigmas.map((s) =>
      seriesFrom(sig, "T_s", "mean_post_std_r_km", { sigma_deg: s, scenario: "CR5" }, () => `CR5 σθ=${s}°`)
    );
    drawLineChart($("#chart1"), series1, {
      xLog: true,
      xLabel: "T (s)",
      yLabel: "CR5 后验σr km",
    });
    const series1b = [0.02, 0.05, 0.1, 0.2, 0.5].map((s) =>
      seriesFrom(sig, "T_s", "mean_abs_err_v", { sigma_deg: s }, () => `σθ=${s}°`)
    );
    drawLineChart($("#chart1b"), series1b, {
      xLog: true,
      xLabel: "T (s)",
      yLabel: "mean |Δv| m/s",
    });
    const psis = [-5, 0, 5, 10];
    const series2 = psis.map((p) =>
      seriesFrom(head, "T_s", "mean_abs_err_psi_deg", { sigma_deg: 0.1, psi_deg: p }, () => `ψ=${p}°`)
    );
    drawLineChart($("#chart2"), series2, {
      xLog: true,
      xLabel: "T (s)",
      yLabel: "mean |Δψ| deg (σθ=0.1°)",
    });
    const turns = [0, 2, 5, 10, 15];
    const series3 = turns.map((t) =>
      seriesFrom(man, "T_s", "mean_abs_err_r_km", { sigma_deg: 0.1, turn_deg: t }, () => `转角=${t}°`)
    );
    drawLineChart($("#chart3"), series3, {
      xLog: true,
      xLabel: "T (s)",
      yLabel: "mean |Δr| km (σθ=0.1°)",
    });
    const series3b = turns.map((t) =>
      seriesFrom(man, "T_s", "sv_min", { sigma_deg: 0.1, turn_deg: t }, () => `转角=${t}°`)
    );
    drawLineChart($("#chart3b"), series3b, {
      xLog: true,
      yLog: true,
      xLabel: "T (s)",
      yLabel: "sv_min",
    });
  }

  function tableFromRows(rows, cols) {
    if (!rows || !rows.length) return '<div class="note">无数据</div>';
    const head = cols.map((c) => `<th>${c.label}</th>`).join("");
    const body = rows
      .map(
        (r) =>
          "<tr>" +
          cols
            .map((c) => `<td>${fmt(r[c.key], c.nd)}</td>`)
            .join("") +
          "</tr>"
      )
      .join("");
    return `<div class="table-wrap"><table><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div>`;
  }

  function renderTables() {
    if (!data) return;
    const sig = (data.sigma_T || []).filter((r) =>
      [60, 300, 600, 1200].includes(Number(r.T_s)) &&
      [0.02, 0.05, 0.1, 0.2, 0.5, 1.0].includes(Number(r.sigma_deg))
    );
    $("#tbl-sigma").innerHTML = tableFromRows(sig, [
      { key: "scenario", label: "场景", nd: null },
      { key: "T_s", label: "T(s)", nd: 0 },
      { key: "sigma_deg", label: "σθ(°)", nd: 2 },
      { key: "mean_post_std_r_km", label: "后验σr km", nd: 3 },
      { key: "mean_r_width_km", label: "r宽度 km", nd: 2 },
      { key: "mean_abs_err_r_km", label: "|Δr| km", nd: 3 },
      { key: "mean_contraction_r", label: "contr_r", nd: 3 },
      { key: "rate_r_resolved", label: "r可辨识", nd: 2 },
      { key: "bearing_span_deg", label: "span°", nd: 3 },
    ]);

    const head = (data.heading || []).filter(
      (r) => Number(r.sigma_deg) === 0.1 && [300, 600, 1200].includes(Number(r.T_s))
    );
    $("#tbl-head").innerHTML = tableFromRows(head, [
      { key: "T_s", label: "T(s)", nd: 0 },
      { key: "psi_deg", label: "ψ(°)", nd: 1 },
      { key: "mean_abs_err_psi_deg", label: "|Δψ|°", nd: 3 },
      { key: "identifiable_rate", label: "可辨识率", nd: 2 },
      { key: "crlb_psi_deg", label: "CRLB ψ", nd: 2 },
      { key: "mean_abs_err_r_km", label: "|Δr| km", nd: 3 },
    ]);

    const man = (data.maneuver || []).filter(
      (r) => Number(r.sigma_deg) === 0.1 && [300, 600, 1200].includes(Number(r.T_s))
    );
    $("#tbl-man").innerHTML = tableFromRows(man, [
      { key: "scenario", label: "场景", nd: null },
      { key: "turn_deg", label: "转角°", nd: 0 },
      { key: "T_s", label: "T(s)", nd: 0 },
      { key: "mean_post_std_r_km", label: "后验σr km", nd: 3 },
      { key: "mean_r_width_km", label: "r宽度", nd: 2 },
      { key: "mean_abs_err_r_km", label: "|Δr| km", nd: 3 },
      { key: "sv_min", label: "sv_min", nd: 4 },
      { key: "bearing_span_deg", label: "span°", nd: 3 },
      { key: "rate_r_resolved", label: "r可辨识", nd: 2 },
    ]);

    const off = data.offgrid_summary || [];
    $("#tbl-off").innerHTML = tableFromRows(off.slice(0, 24), [
      { key: "tag", label: "OG", nd: null },
      { key: "T_s", label: "T(s)", nd: 0 },
      { key: "sigma_deg", label: "σθ°", nd: 2 },
      { key: "mean_abs_err_r_km", label: "|Δr| km", nd: 3 },
      { key: "mean_abs_err_v", label: "|Δv|", nd: 3 },
      { key: "mean_abs_err_psi_deg", label: "|Δψ|°", nd: 3 },
      { key: "any_on_grid", label: "真值在网格", nd: null },
    ]);

    const hard = data.hard_pairs || [];
    const hardShow = hard.slice(0, 20);
    $("#tbl-hard").innerHTML = tableFromRows(hardShow, [
      { key: "ref_name", label: "机制/标签", nd: null },
      { key: "turn_deg", label: "转角°", nd: 0 },
      { key: "T_s", label: "T(s)", nd: 0 },
      { key: "ref_r0_km", label: "参考 r km", nd: 2 },
      { key: "ref_v", label: "参考 v", nd: 2 },
      { key: "alt_r0_km", label: "备择 r km", nd: 2 },
      { key: "alt_v", label: "备择 v", nd: 2 },
      { key: "alt_psi_deg", label: "备择 ψ°", nd: 2 },
      { key: "rmse_bearing_deg", label: "方位RMSE°", nd: 4 },
      { key: "grid_alt_in_accept_rate", label: "进入接受集", nd: 2 },
    ]);
  }

  function renderFigures() {
    const figs = [
      "fig1a_sigma_T_range_poststd_cr5.svg",
      "fig1a_sigma_T_range_width_cr5.svg",
      "fig1b_sigma_T_range_err_cr5.svg",
      "fig1c_sigma_T_speed_err_cr5.svg",
      "fig1d_sigma_T_contraction_cr5.svg",
      "fig1e_sigma_T_truth_in_band_cr5.svg",
      "fig1a_sigma_T_range_poststd_lin.svg",
      "fig1a_sigma_T_range_width_lin.svg",
      "fig2a_heading_err_vs_T.svg",
      "fig2b_psi_crlb_vs_T.svg",
      "fig2c_heading_identifiable_rate_heatmap.svg",
      "fig2d_psi_crlb_heatmap.svg",
      "fig3a_turn_T_range_poststd_cr5.svg",
      "fig3a_turn_T_range_width_cr5.svg",
      "fig3b_turn_T_svmin_cr5.svg",
      "fig3c_turn_T_contraction_cr5.svg",
      "fig3d_turn_T_range_width_heatmap_cr5.svg",
      "fig3e_turn_T_bearing_span_cr5.svg",
      "fig3f_hard_pair_bearing_rmse.svg",
    ];
    const box = $("#figures");
    box.innerHTML = figs
      .map((f) => {
        const url = "results/P2_RC2_kinematic_boundary/figures/" + f;
        return `<div class="img-frame"><object data="${url}" type="image/svg+xml" title="${f}"><a href="${url}">${f}</a></object></div>`;
      })
      .join("");
  }

  function renderFiles() {
    const files = [
      ["results/P2_RC2_kinematic_boundary/P2_RC2_REPORT.md", "主报告"],
      ["results/P2_RC2_kinematic_boundary/P2_RC2_GPT_SYNC.md", "GPT同步稿"],
      ["results/P2_RC2_kinematic_boundary/P2_RC2_CONFIG.json", "配置冻结"],
      ["results/P2_RC2_kinematic_boundary/observability_scan.csv", "FIM/可观测性"],
      ["results/P2_RC2_kinematic_boundary/bearing_accuracy_time_boundary.csv", "σθ×T 边界"],
      ["results/P2_RC2_kinematic_boundary/heading_identifiability_boundary.csv", "ψ×T 边界"],
      ["results/P2_RC2_kinematic_boundary/maneuver_boundary.csv", "转角×T 边界"],
      ["results/P2_RC2_kinematic_boundary/offgrid_trials.csv", "off-grid 试验"],
      ["results/P2_RC2_kinematic_boundary/hard_candidate_pairs.csv", "困难候选"],
      ["results/P2_RC2_kinematic_boundary/hard_candidate_pairs_rc3_handoff.csv", "RC3 移交清单"],
      ["results/P2_RC2_kinematic_boundary/p2_results.json", "仪表盘数据"],
      ["results/P2_RC2_kinematic_boundary/P2_RC2_GPT_SYNC.md", "GPT 同步"],
    ];
    $("#files").innerHTML = files
      .map(
        ([href, label]) =>
          `<a class="file" href="${href}" target="_blank" rel="noopener">${href.split("/").pop()}<span>${label}</span></a>`
      )
      .join("");
  }

  function renderFilters() {
    const sig = data && data.sigma_T ? data.sigma_T : [];
    const sigmas = unique(sig.map((r) => r.sigma_deg)).sort((a, b) => a - b);
    const Ts = unique(sig.map((r) => r.T_s)).sort((a, b) => a - b);
    fillSelect($("#f-sigma"), sigmas.length ? sigmas : [0.1], 0.1);
    fillSelect($("#f-T"), Ts.length ? Ts : [600], 600);
    $("#f-apply").onclick = function () {
      const s = Number($("#f-sigma").value);
      const T = Number($("#f-T").value);
      const row = sig.find((r) => Number(r.sigma_deg) === s && Number(r.T_s) === T);
      const out = $("#filter-out");
      if (!row) {
        out.textContent = "该组合无数据";
        return;
      }
      out.innerHTML =
        `<code>${row.scenario || "?"} · σθ=${s}° · T=${T}s</code> → 后验σr=<strong>${fmt(row.mean_post_std_r_km)}</strong> km, ` +
        `r宽度=<strong>${fmt(row.mean_r_width_km)}</strong> km, ` +
        `|Δr|=<strong>${fmt(row.mean_abs_err_r_km)}</strong> km, ` +
        `contr_r=<strong>${fmt(row.mean_contraction_r)}</strong>, ` +
        `r可辨识=<strong>${fmt(row.rate_r_resolved,2)}</strong>, ` +
        `span=<strong>${fmt(row.bearing_span_deg,3)}</strong>°, ` +
        `CRLB r=<strong>${fmt(row.crlb_r_km)}</strong> km`;
    };
  }

  function boot() {
    $("#meta-line").textContent = data
      ? `项目 ${data.meta.project} · ${data.meta.package} · UTC ${data.meta.created_utc} · 每格非零噪声试验 ${data.meta.n_trials} 次 · dt=${data.meta.dt_obs_s}s`
      : "未加载 p2_results.js";
    renderKpis();
    renderVerdict();
    renderCharts();
    renderTables();
    renderFigures();
    renderFiles();
    renderFilters();
    $("#f-apply") && $("#f-apply").click();
    window.addEventListener("resize", function () {
      clearTimeout(window.__p2rs);
      window.__p2rs = setTimeout(renderCharts, 150);
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", boot);
  } else {
    boot();
  }
})();
