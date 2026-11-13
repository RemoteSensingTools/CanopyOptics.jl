/* Interactive D3.js v7 figure engine for CanopyOptics.jl documentation.
 * Supports multi-panel Cartesian charts (lines, bands, grouped bars, scatter,
 * dual y-axes, reference lines/spans, insets, clickable legend toggles),
 * 2D parameter-slider heatmaps, and 2D polar hemisphere BRDF/NDVI surfaces.
 */
(function () {
  "use strict";

  const PALETTE = [
    "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e",
    "#9467bd", "#0072B2", "#009E73", "#CC79A7"
  ];

  const DASH_MAP = {
    solid: null,
    dash: "6,4",
    dot: "2,3",
    dashdot: "7,3,2,3"
  };

  let sharedTooltip = null;
  function getTooltip() {
    if (!sharedTooltip) {
      sharedTooltip = document.createElement("div");
      sharedTooltip.className = "canopy-d3-tooltip";
      document.body.appendChild(sharedTooltip);
    }
    return sharedTooltip;
  }

  function showTooltip(evt, html) {
    const tt = getTooltip();
    tt.innerHTML = html;
    tt.classList.add("is-visible");
    const pad = 14;
    const vw = window.innerWidth;
    const vh = window.innerHeight;
    const rect = tt.getBoundingClientRect();
    let left = evt.clientX + pad;
    let top = evt.clientY + pad;
    if (left + rect.width + 10 > vw) left = Math.max(8, evt.clientX - rect.width - pad);
    if (top + rect.height + 10 > vh) top = Math.max(8, evt.clientY - rect.height - pad);
    tt.style.left = left + "px";
    tt.style.top = top + "px";
  }

  function hideTooltip() {
    if (sharedTooltip) sharedTooltip.classList.remove("is-visible");
  }

  function fmtNum(v) {
    if (v === null || v === undefined || Number.isNaN(v)) return "—";
    const av = Math.abs(v);
    if (av === 0) return "0";
    if (av >= 1000 || av < 1e-3) return v.toExponential(2);
    if (av >= 100) return v.toFixed(1);
    if (av >= 10) return v.toFixed(2);
    if (av >= 1) return v.toFixed(3);
    return v.toFixed(4);
  }

  function getColormap(name) {
    if (name === "balance" || name === "rdbu") {
      return function (t) {
        const u = Math.max(0, Math.min(1, t));
        if (window.d3 && d3.interpolateRdBu) return d3.interpolateRdBu(1 - u);
        const r = Math.round(255 * Math.min(1, 2 * u));
        const b = Math.round(255 * Math.min(1, 2 * (1 - u)));
        const g = Math.round(245 * (1 - Math.abs(2 * u - 1)));
        return "rgb(" + r + "," + g + "," + b + ")";
      };
    }
    return function (t) {
      const u = Math.max(0, Math.min(1, t));
      if (window.d3 && d3.interpolateViridis) return d3.interpolateViridis(u);
      return "rgb(" + Math.round(68 + 185 * u) + "," + Math.round(1 + 220 * u) + "," + Math.round(84 + 80 * (1 - u)) + ")";
    };
  }

  function parseRGB(cssColor) {
    if (window.d3 && d3.color) {
      const c = d3.color(cssColor);
      if (c) return [c.r | 0, c.g | 0, c.b | 0];
    }
    return [100, 120, 180];
  }

  function buildColorLUT(name) {
    const fn = getColormap(name);
    const lut = new Uint8Array(256 * 3);
    for (let i = 0; i < 256; i++) {
      const rgb = parseRGB(fn(i / 255));
      lut[3 * i] = rgb[0];
      lut[3 * i + 1] = rgb[1];
      lut[3 * i + 2] = rgb[2];
    }
    return lut;
  }

  function drawSymbol(sel, kind, size) {
    const r = size || 3.4;
    if (kind === "square") {
      sel.append("rect")
        .attr("x", -r).attr("y", -r)
        .attr("width", 2 * r).attr("height", 2 * r);
    } else if (kind === "triangle-up") {
      const h = r * 1.25;
      sel.append("polygon").attr("points", `0,${-h} ${-h},${h} ${h},${h}`);
    } else if (kind === "triangle-down") {
      const h = r * 1.25;
      sel.append("polygon").attr("points", `0,${h} ${-h},${-h} ${h},${-h}`);
    } else if (kind === "cross") {
      const h = r * 1.1;
      sel.append("path")
        .attr("d", `M${-h},${-h}L${h},${h}M${-h},${h}L${h},${-h}`)
        .attr("stroke-width", 1.6)
        .attr("fill", "none");
    } else {
      sel.append("circle").attr("r", r);
    }
  }

  function starPoints(cx, cy, outerR, innerR) {
    const pts = [];
    for (let i = 0; i < 10; i++) {
      const ang = -Math.PI / 2 + (i * Math.PI) / 5;
      const rad = i % 2 === 0 ? outerR : innerR;
      pts.push((cx + rad * Math.cos(ang)).toFixed(2) + "," + (cy + rad * Math.sin(ang)).toFixed(2));
    }
    return pts.join(" ");
  }

  function styleAxes(svg) {
    svg.selectAll(".cp-axis path, .cp-axis line")
      .attr("stroke", "var(--cp-axis)")
      .attr("stroke-opacity", 0.65);
    svg.selectAll(".cp-axis text")
      .attr("fill", "var(--cp-text)")
      .style("font-size", "11px");
    svg.selectAll(".cp-grid line")
      .attr("stroke", "var(--cp-grid)")
      .attr("stroke-dasharray", "2,2");
    svg.selectAll(".cp-grid path").attr("stroke", "none");
  }

  function renderCartesianPanel(panelEl, panel, rowHeight) {
    const series = panel.series || [];
    const muted = new Set();

    const legendSeries = series.filter((s) => s.label && s.showInLegend !== false);
    let legendEl = null;
    if (legendSeries.length > 0) {
      legendEl = document.createElement("div");
      legendEl.className = "canopy-d3-legend";
      panelEl.appendChild(legendEl);
    }

    const stage = document.createElement("div");
    stage.className = "canopy-d3-stage";
    panelEl.appendChild(stage);

    function draw() {
      stage.innerHTML = "";
      const W = Math.max(260, stage.clientWidth || 440);
      const H = panel.height || rowHeight || 300;
      const hasY2 = Boolean(panel.y2label) || series.some((s) => s.yAxis === 2);
      const cats = panel.categories || null;
      const maxCatLines = cats
        ? Math.max(...cats.map((c) => String(c).split("\n").length))
        : 1;
      const margin = {
        top: 12,
        right: hasY2 ? 54 : 18,
        bottom: (panel.xlabel ? 40 : 26) + (maxCatLines - 1) * 12,
        left: 56
      };
      const iw = Math.max(80, W - margin.left - margin.right);
      const ih = Math.max(80, H - margin.top - margin.bottom);

      const svg = d3.select(stage)
        .append("svg")
        .attr("viewBox", `0 0 ${W} ${H}`)
        .attr("preserveAspectRatio", "xMidYMid meet");

      const g = svg.append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);

      const activeSeries = series.filter((_, idx) => !muted.has(idx));

      // X scale
      let xScale;
      if (cats) {
        xScale = d3.scaleBand()
          .domain(cats.map((_, i) => i))
          .range([0, iw])
          .paddingInner(0.22)
          .paddingOuter(0.12);
      } else {
        let xmin = Infinity, xmax = -Infinity;
        if (panel.xlim) {
          xmin = panel.xlim[0];
          xmax = panel.xlim[1];
        } else {
          for (const s of activeSeries.length ? activeSeries : series) {
            if (s.x) {
              for (const xv of s.x) {
                if (xv !== null && Number.isFinite(xv)) {
                  if (xv < xmin) xmin = xv;
                  if (xv > xmax) xmax = xv;
                }
              }
            }
          }
          if (!Number.isFinite(xmin)) { xmin = 0; xmax = 1; }
          if (xmin === xmax) { xmin -= 0.5; xmax += 0.5; }
        }
        xScale = d3.scaleLinear().domain([xmin, xmax]).range([0, iw]);
      }

      // Y1 & Y2 scales
      function computeYDomain(axisIdx, lim) {
        if (lim) return lim;
        let ymin = Infinity, ymax = -Infinity;
        const pool = (activeSeries.length ? activeSeries : series).filter(
          (s) => (s.yAxis || 1) === axisIdx
        );
        for (const s of pool) {
          if (s.kind === "bar") ymin = Math.min(ymin, 0);
          const arrs = s.kind === "band" ? [s.y0 || [], s.y1 || []] : [s.y || []];
          for (const arr of arrs) {
            for (const yv of arr) {
              if (yv !== null && Number.isFinite(yv)) {
                if (yv < ymin) ymin = yv;
                if (yv > ymax) ymax = yv;
              }
            }
          }
        }
        if (axisIdx === 1 && panel.hlines) {
          for (const hl of panel.hlines) {
            if (Number.isFinite(hl.y)) {
              if (hl.y < ymin) ymin = hl.y;
              if (hl.y > ymax) ymax = hl.y;
            }
          }
        }
        if (!Number.isFinite(ymin)) return [0, 1];
        if (ymin === ymax) return [ymin - 0.5, ymax + 0.5];
        const span = ymax - ymin;
        const lo = sHasOnlyNonNegBars(pool) && ymin >= 0 ? 0 : ymin - 0.06 * span;
        return [lo, ymax + 0.08 * span];
      }

      function sHasOnlyNonNegBars(pool) {
        return pool.length > 0 && pool.every((s) => s.kind === "bar");
      }

      const y1Scale = d3.scaleLinear().domain(computeYDomain(1, panel.ylim)).range([ih, 0]);
      if (!panel.ylim) y1Scale.nice(5);
      const y2Scale = hasY2
        ? d3.scaleLinear().domain(computeYDomain(2, panel.y2lim)).range([ih, 0])
        : null;
      if (y2Scale && !panel.y2lim) y2Scale.nice(5);

      // Gridlines
      g.append("g")
        .attr("class", "cp-grid")
        .call(d3.axisLeft(y1Scale).ticks(5).tickSize(-iw).tickFormat(""));

      if (!cats) {
        const xGridAxis = d3.axisBottom(xScale).ticks(panel.xticks ? panel.xticks.length : 6).tickSize(-ih).tickFormat("");
        if (panel.xticks) xGridAxis.tickValues(panel.xticks);
        g.append("g")
          .attr("class", "cp-grid")
          .attr("transform", `translate(0,${ih})`)
          .call(xGridAxis);
      }

      // Clip path
      const clipId = "cp-clip-" + Math.random().toString(36).slice(2, 9);
      svg.append("defs").append("clipPath")
        .attr("id", clipId)
        .append("rect")
        .attr("x", 0).attr("y", 0)
        .attr("width", iw).attr("height", ih);

      const plotArea = g.append("g").attr("clip-path", `url(#${clipId})`);

      // Vertical spans
      if (panel.vspans && !cats) {
        for (const sp of panel.vspans) {
          const x0 = xScale(sp.x0);
          const x1 = xScale(sp.x1);
          plotArea.append("rect")
            .attr("x", Math.min(x0, x1))
            .attr("y", 0)
            .attr("width", Math.abs(x1 - x0))
            .attr("height", ih)
            .attr("fill", sp.color || "#F0E442")
            .attr("fill-opacity", sp.opacity !== undefined ? sp.opacity : 0.18);
        }
      }

      // Reference hlines & vlines
      if (panel.hlines) {
        for (const hl of panel.hlines) {
          const yv = y1Scale(hl.y);
          plotArea.append("line")
            .attr("x1", 0).attr("x2", iw)
            .attr("y1", yv).attr("y2", yv)
            .attr("stroke", hl.color || "var(--cp-axis)")
            .attr("stroke-width", hl.width || 1.1)
            .attr("stroke-dasharray", DASH_MAP[hl.dash || "dash"] || null);
        }
      }
      if (panel.vlines) {
        for (const vl of panel.vlines) {
          const xv = cats
            ? (xScale(Math.floor(vl.x)) || 0) + xScale.bandwidth() * (vl.x - Math.floor(vl.x) + 0.6)
            : xScale(vl.x);
          plotArea.append("line")
            .attr("x1", xv).attr("x2", xv)
            .attr("y1", 0).attr("y2", ih)
            .attr("stroke", vl.color || "var(--cp-axis)")
            .attr("stroke-width", vl.width || 1.1)
            .attr("stroke-dasharray", DASH_MAP[vl.dash || "dash"] || null);
        }
      }

      // Grouped bars
      const barEntries = series
        .map((s, idx) => ({ s, idx }))
        .filter((e) => e.s.kind === "bar" && !muted.has(e.idx));
      if (cats && barEntries.length > 0) {
        const subScale = d3.scaleBand()
          .domain(d3.range(barEntries.length))
          .range([0, xScale.bandwidth()])
          .padding(0.08);
        barEntries.forEach(({ s, idx }, bIdx) => {
          const ys = s.yAxis === 2 && y2Scale ? y2Scale : y1Scale;
          const yZero = ys(0);
          const col = s.color || PALETTE[idx % PALETTE.length];
          cats.forEach((catName, cIdx) => {
            const val = s.y[cIdx];
            if (val === null || val === undefined || !Number.isFinite(val)) return;
            const yTop = Math.min(yZero, ys(val));
            const bh = Math.max(0, Math.abs(ys(val) - yZero));
            plotArea.append("rect")
              .attr("x", xScale(cIdx) + subScale(bIdx))
              .attr("y", yTop)
              .attr("width", subScale.bandwidth())
              .attr("height", bh)
              .attr("fill", col)
              .attr("fill-opacity", s.opacity !== undefined ? s.opacity : 0.88)
              .attr("rx", 1.5)
              .style("pointer-events", "all")
              .on("mousemove", (evt) => {
                const cleanCat = String(catName).replace(/\n/g, " ");
                showTooltip(
                  evt,
                  `<div class="canopy-d3-tooltip-head">${cleanCat}</div>` +
                  `<div class="canopy-d3-tooltip-row"><span><span class="canopy-d3-tooltip-dot" style="background:${col}"></span>${s.label || "Series"}</span><strong>${fmtNum(val)}</strong></div>`
                );
              })
              .on("mouseleave", hideTooltip);
          });
        });
      }

      // Bands, Lines, Scatters
      series.forEach((s, idx) => {
        if (muted.has(idx) || s.kind === "bar") return;
        const col = s.color || PALETTE[idx % PALETTE.length];
        const ys = s.yAxis === 2 && y2Scale ? y2Scale : y1Scale;

        if (s.kind === "band") {
          const area = d3.area()
            .defined((_, i) => Number.isFinite(s.x[i]) && Number.isFinite(s.y0[i]) && Number.isFinite(s.y1[i]))
            .x((_, i) => xScale(s.x[i]))
            .y0((_, i) => ys(s.y0[i]))
            .y1((_, i) => ys(s.y1[i]));
          plotArea.append("path")
            .datum(d3.range(s.x.length))
            .attr("d", area)
            .attr("fill", col)
            .attr("fill-opacity", s.opacity !== undefined ? s.opacity : 0.18);
          return;
        }

        if (s.kind === "scatter") {
          const pts = plotArea.append("g")
            .attr("fill", col)
            .attr("stroke", col)
            .attr("fill-opacity", s.opacity !== undefined ? s.opacity : 0.68);
          for (let i = 0; i < s.x.length; i++) {
            if (!Number.isFinite(s.x[i]) || !Number.isFinite(s.y[i])) continue;
            const pg = pts.append("g")
              .attr("transform", `translate(${xScale(s.x[i])},${ys(s.y[i])})`)
              .style("pointer-events", "all")
              .on("mousemove", (evt) => {
                showTooltip(
                  evt,
                  `<div class="canopy-d3-tooltip-head">${s.label || "Point"}</div>` +
                  `<div class="canopy-d3-tooltip-row"><span>x</span><strong>${fmtNum(s.x[i])}</strong></div>` +
                  `<div class="canopy-d3-tooltip-row"><span>y</span><strong>${fmtNum(s.y[i])}</strong></div>`
                );
              })
              .on("mouseleave", hideTooltip);
            drawSymbol(pg, s.marker || "circle", s.size || 3.2);
          }
          return;
        }

        // Default: line
        if (s.fillZero) {
          const area0 = d3.area()
            .defined((_, i) => Number.isFinite(s.x[i]) && Number.isFinite(s.y[i]))
            .x((_, i) => xScale(s.x[i]))
            .y0(ys(0))
            .y1((_, i) => ys(s.y[i]));
          plotArea.append("path")
            .datum(d3.range(s.x.length))
            .attr("d", area0)
            .attr("fill", col)
            .attr("fill-opacity", s.fillOpacity !== undefined ? s.fillOpacity : 0.14);
        }

        const lineGen = d3.line()
          .defined((_, i) => Number.isFinite(s.x[i]) && Number.isFinite(s.y[i]))
          .x((_, i) => xScale(s.x[i]))
          .y((_, i) => ys(s.y[i]));

        plotArea.append("path")
          .datum(d3.range(s.x.length))
          .attr("d", lineGen)
          .attr("fill", "none")
          .attr("stroke", col)
          .attr("stroke-width", s.width || 2.0)
          .attr("stroke-opacity", s.opacity !== undefined ? s.opacity : 1.0)
          .attr("stroke-dasharray", DASH_MAP[s.dash || "solid"] || null);

        if (s.marker) {
          const mg = plotArea.append("g")
            .attr("fill", col)
            .attr("stroke", col)
            .attr("opacity", s.opacity !== undefined ? s.opacity : 1.0);
          for (let i = 0; i < s.x.length; i++) {
            if (!Number.isFinite(s.x[i]) || !Number.isFinite(s.y[i])) continue;
            const pg = mg.append("g")
              .attr("transform", `translate(${xScale(s.x[i])},${ys(s.y[i])})`);
            drawSymbol(pg, s.marker, s.markerSize || 3.3);
          }
        }
      });

      // Axes
      const xAxis = d3.axisBottom(xScale);
      if (cats) {
        xAxis.tickFormat((i) => cats[i] || "");
      } else if (panel.xticks) {
        xAxis.tickValues(panel.xticks);
      } else {
        xAxis.ticks(6);
      }
      const gx = g.append("g")
        .attr("class", "cp-axis")
        .attr("transform", `translate(0,${ih})`)
        .call(xAxis);

      if (cats) {
        gx.selectAll(".tick text").each(function (i) {
          const el = d3.select(this);
          const lines = String(cats[i] || "").split("\n");
          if (lines.length > 1) {
            el.text(null);
            lines.forEach((ln, li) => {
              el.append("tspan")
                .attr("x", 0)
                .attr("dy", li === 0 ? "0.75em" : "1.12em")
                .text(ln);
            });
          }
        }).style("font-size", "10px");
      }

      g.append("g")
        .attr("class", "cp-axis")
        .call(d3.axisLeft(y1Scale).ticks(5));

      if (hasY2 && y2Scale) {
        const gy2 = g.append("g")
          .attr("class", "cp-axis")
          .attr("transform", `translate(${iw},0)`)
          .call(d3.axisRight(y2Scale).ticks(5));
        if (panel.y2color) {
          gy2.selectAll("text").style("fill", panel.y2color);
        }
      }

      styleAxes(svg);

      // Axis labels
      if (panel.xlabel) {
        g.append("text")
          .attr("x", iw / 2)
          .attr("y", ih + margin.bottom - 6)
          .attr("text-anchor", "middle")
          .attr("fill", "var(--cp-text)")
          .style("font-size", "11.5px")
          .text(panel.xlabel);
      }
      if (panel.ylabel) {
        g.append("text")
          .attr("transform", "rotate(-90)")
          .attr("x", -ih / 2)
          .attr("y", -margin.left + 15)
          .attr("text-anchor", "middle")
          .attr("fill", panel.y1color || "var(--cp-text)")
          .style("font-size", "11.5px")
          .text(panel.ylabel);
      }
      if (panel.y2label && hasY2) {
        g.append("text")
          .attr("transform", "rotate(90)")
          .attr("x", ih / 2)
          .attr("y", -iw - margin.right + 14)
          .attr("text-anchor", "middle")
          .attr("fill", panel.y2color || "var(--cp-text)")
          .style("font-size", "11.5px")
          .text(panel.y2label);
      }

      // Optional inset (e.g. Fig 7A refractive index dispersion)
      if (panel.inset) {
        const ins = panel.inset;
        const [fx, fy, fw, fh] = ins.pos || [0.55, 0.06, 0.41, 0.36];
        const ix0 = fx * iw, iy0 = fy * ih, iiw = fw * iw, iih = fh * ih;
        const ig = g.append("g").attr("transform", `translate(${ix0},${iy0})`);
        ig.append("rect")
          .attr("x", -26).attr("y", -14)
          .attr("width", iiw + 32).attr("height", iih + 30)
          .attr("rx", 5)
          .attr("fill", "var(--cp-bg)")
          .attr("fill-opacity", 0.92)
          .attr("stroke", "var(--cp-border)");
        const ixDom = ins.xlim || [400, 2450];
        const ixScale = d3.scaleLinear().domain(ixDom).range([0, iiw]);
        const iyScale = d3.scaleLinear().domain(ins.ylim || [1.25, 1.50]).range([iih, 0]);
        if (ins.title) {
          ig.append("text")
            .attr("x", iiw / 2).attr("y", -3)
            .attr("text-anchor", "middle")
            .attr("fill", "var(--cp-text)")
            .style("font-size", "9.5px")
            .style("font-weight", "600")
            .text(ins.title);
        }
        (ins.series || []).forEach((s, k) => {
          const col = s.color || PALETTE[k % PALETTE.length];
          const lg = d3.line()
            .defined((_, i) => Number.isFinite(s.x[i]) && Number.isFinite(s.y[i]) && s.x[i] >= ixDom[0] && s.x[i] <= ixDom[1])
            .x((_, i) => ixScale(s.x[i]))
            .y((_, i) => iyScale(s.y[i]));
          ig.append("path")
            .datum(d3.range(s.x.length))
            .attr("d", lg)
            .attr("fill", "none")
            .attr("stroke", col)
            .attr("stroke-width", s.width || 1.5)
            .attr("stroke-dasharray", DASH_MAP[s.dash || "solid"] || null);
        });
        ig.append("g")
          .attr("class", "cp-axis")
          .attr("transform", `translate(0,${iih})`)
          .call(d3.axisBottom(ixScale).ticks(3));
        ig.append("g")
          .attr("class", "cp-axis")
          .call(d3.axisLeft(iyScale).ticks(3));
        ig.selectAll(".cp-axis text").style("font-size", "8.5px").attr("fill", "var(--cp-text)");
        ig.selectAll(".cp-axis path, .cp-axis line").attr("stroke", "var(--cp-axis)");
      }

      // Hover crosshair for continuous line charts
      const lineEntries = series
        .map((s, idx) => ({ s, idx }))
        .filter((e) => (!e.s.kind || e.s.kind === "line") && !muted.has(e.idx) && e.s.x && e.s.x.length > 0);
      if (!cats && lineEntries.length > 0) {
        const hoverLine = g.append("line")
          .attr("y1", 0).attr("y2", ih)
          .attr("stroke", "var(--cp-axis)")
          .attr("stroke-width", 1)
          .attr("stroke-dasharray", "3,3")
          .style("opacity", 0)
          .style("pointer-events", "none");

        const hoverDots = g.append("g").style("pointer-events", "none");

        const bisect = d3.bisector((d) => d).center;
        g.append("rect")
          .attr("x", 0).attr("y", 0)
          .attr("width", iw).attr("height", ih)
          .attr("fill", "transparent")
          .style("cursor", "crosshair")
          .on("mousemove", function (evt) {
            const [mx] = d3.pointer(evt, this);
            const xVal = xScale.invert(mx);
            hoverLine.attr("x1", mx).attr("x2", mx).style("opacity", 0.75);
            hoverDots.selectAll("*").remove();

            let rowsHtml = "";
            let refX = xVal;
            lineEntries.forEach(({ s, idx }, k) => {
              const i = bisect(s.x, xVal);
              if (i < 0 || i >= s.x.length || !Number.isFinite(s.y[i])) return;
              if (k === 0) refX = s.x[i];
              const col = s.color || PALETTE[idx % PALETTE.length];
              const ys = s.yAxis === 2 && y2Scale ? y2Scale : y1Scale;
              hoverDots.append("circle")
                .attr("cx", xScale(s.x[i]))
                .attr("cy", ys(s.y[i]))
                .attr("r", 3.6)
                .attr("fill", col)
                .attr("stroke", "var(--cp-bg)")
                .attr("stroke-width", 1.2);
              if (s.label) {
                rowsHtml += `<div class="canopy-d3-tooltip-row"><span><span class="canopy-d3-tooltip-dot" style="background:${col}"></span>${s.label}</span><strong>${fmtNum(s.y[i])}</strong></div>`;
              }
            });
            if (rowsHtml) {
              const xName = panel.xlabel ? panel.xlabel.split("(")[0].trim() : "x";
              showTooltip(evt, `<div class="canopy-d3-tooltip-head">${xName}: ${fmtNum(refX)}</div>` + rowsHtml);
            }
          })
          .on("mouseleave", () => {
            hoverLine.style("opacity", 0);
            hoverDots.selectAll("*").remove();
            hideTooltip();
          });
      }
    }

    if (legendEl) {
      series.forEach((s, idx) => {
        if (!s.label || s.showInLegend === false) return;
        const item = document.createElement("span");
        item.className = "canopy-d3-legend-item";
        item.title = "Click to toggle series";
        const col = s.color || PALETTE[idx % PALETTE.length];
        const sw = document.createElementNS("http://www.w3.org/2000/svg", "svg");
        sw.setAttribute("class", "canopy-d3-swatch");
        sw.setAttribute("viewBox", "0 0 18 10");
        if (s.kind === "bar" || s.kind === "band") {
          const r = document.createElementNS("http://www.w3.org/2000/svg", "rect");
          r.setAttribute("x", "2"); r.setAttribute("y", "1.5");
          r.setAttribute("width", "14"); r.setAttribute("height", "7");
          r.setAttribute("rx", "1.5");
          r.setAttribute("fill", col);
          sw.appendChild(r);
        } else if (s.kind === "scatter") {
          const c = document.createElementNS("http://www.w3.org/2000/svg", "circle");
          c.setAttribute("cx", "9"); c.setAttribute("cy", "5"); c.setAttribute("r", "3.2");
          c.setAttribute("fill", col);
          sw.appendChild(c);
        } else {
          const ln = document.createElementNS("http://www.w3.org/2000/svg", "line");
          ln.setAttribute("x1", "0"); ln.setAttribute("x2", "18");
          ln.setAttribute("y1", "5"); ln.setAttribute("y2", "5");
          ln.setAttribute("stroke", col);
          ln.setAttribute("stroke-width", "2.2");
          const da = DASH_MAP[s.dash || "solid"];
          if (da) ln.setAttribute("stroke-dasharray", da);
          sw.appendChild(ln);
        }
        const txt = document.createElement("span");
        txt.textContent = s.label;
        item.appendChild(sw);
        item.appendChild(txt);
        item.addEventListener("click", () => {
          if (muted.has(idx)) {
            muted.delete(idx);
            item.classList.remove("is-muted");
          } else {
            muted.add(idx);
            item.classList.add("is-muted");
          }
          draw();
        });
        legendEl.appendChild(item);
      });
    }

    draw();
    return draw;
  }

  function appendColorbar(svg, x0, y0, cbW, cbH, zmin, zmax, colormapName, label) {
    const gradId = "cp-cb-" + Math.random().toString(36).slice(2, 9);
    const defs = svg.append("defs");
    const grad = defs.append("linearGradient")
      .attr("id", gradId)
      .attr("x1", "0%").attr("y1", "100%")
      .attr("x2", "0%").attr("y2", "0%");
    const cmap = getColormap(colormapName);
    for (let k = 0; k <= 16; k++) {
      const t = k / 16;
      grad.append("stop")
        .attr("offset", (t * 100).toFixed(1) + "%")
        .attr("stop-color", cmap(t));
    }
    const cg = svg.append("g").attr("transform", `translate(${x0},${y0})`);
    cg.append("rect")
      .attr("x", 0).attr("y", 0)
      .attr("width", cbW).attr("height", cbH)
      .attr("fill", `url(#${gradId})`)
      .attr("stroke", "var(--cp-border)");
    const cbScale = d3.scaleLinear().domain([zmin, zmax]).range([cbH, 0]);
    cg.append("g")
      .attr("class", "cp-axis")
      .attr("transform", `translate(${cbW},0)`)
      .call(d3.axisRight(cbScale).ticks(5));
    if (label) {
      cg.append("text")
        .attr("x", cbW / 2)
        .attr("y", -8)
        .attr("text-anchor", "middle")
        .attr("fill", "var(--cp-text)")
        .style("font-size", "10.5px")
        .style("font-weight", "600")
        .text(label);
    }
  }

  function renderHeatmapPanel(panelEl, panel, rowHeight) {
    const stage = document.createElement("div");
    stage.className = "canopy-d3-stage";
    panelEl.appendChild(stage);

    function draw() {
      stage.innerHTML = "";
      const W = Math.max(260, stage.clientWidth || 440);
      const H = panel.height || rowHeight || 300;
      const margin = { top: 20, right: 78, bottom: 42, left: 54 };
      const iw = Math.max(80, W - margin.left - margin.right);
      const ih = Math.max(80, H - margin.top - margin.bottom);

      const xArr = panel.x;
      const yArr = panel.y;
      const zMat = panel.z; // [ny][nx]
      const ny = yArr.length;
      const nx = xArr.length;

      let zmin = Infinity, zmax = -Infinity;
      if (panel.clim) {
        zmin = panel.clim[0];
        zmax = panel.clim[1];
      } else {
        for (let j = 0; j < ny; j++) {
          for (let i = 0; i < nx; i++) {
            const v = zMat[j][i];
            if (Number.isFinite(v)) {
              if (v < zmin) zmin = v;
              if (v > zmax) zmax = v;
            }
          }
        }
        if (!Number.isFinite(zmin)) { zmin = 0; zmax = 1; }
        if (zmin === zmax) { zmax = zmin + 1e-6; }
      }

      // Offscreen raster -> dataURL embedded inside SVG image so viewBox scaling stays 100% aligned!
      const lut = buildColorLUT(panel.colormap || "viridis");
      const cw = Math.max(nx * 4, 160);
      const ch = Math.max(ny * 4, 160);
      const cvs = document.createElement("canvas");
      cvs.width = cw;
      cvs.height = ch;
      const ctx = cvs.getContext("2d");
      const img = ctx.createImageData(cw, ch);
      const xmin = xArr[0], xmax = xArr[nx - 1];
      const ymin = yArr[0], ymax = yArr[ny - 1];
      const bisectX = d3.bisector((d) => d).right;
      const bisectY = d3.bisector((d) => d).right;

      for (let py = 0; py < ch; py++) {
        const yv = ymax - ((py + 0.5) / ch) * (ymax - ymin);
        let j1 = bisectY(yArr, yv);
        j1 = Math.max(1, Math.min(ny - 1, j1));
        const j0 = j1 - 1;
        const ty = yArr[j1] === yArr[j0] ? 0 : Math.max(0, Math.min(1, (yv - yArr[j0]) / (yArr[j1] - yArr[j0])));
        for (let px = 0; px < cw; px++) {
          const xv = xmin + ((px + 0.5) / cw) * (xmax - xmin);
          let i1 = bisectX(xArr, xv);
          i1 = Math.max(1, Math.min(nx - 1, i1));
          const i0 = i1 - 1;
          const tx = xArr[i1] === xArr[i0] ? 0 : Math.max(0, Math.min(1, (xv - xArr[i0]) / (xArr[i1] - xArr[i0])));
          const v00 = zMat[j0][i0], v10 = zMat[j0][i1];
          const v01 = zMat[j1][i0], v11 = zMat[j1][i1];
          const val = (1 - ty) * ((1 - tx) * v00 + tx * v10) + ty * ((1 - tx) * v01 + tx * v11);
          const u = Math.max(0, Math.min(255, Math.round(255 * (val - zmin) / (zmax - zmin))));
          const off = 4 * (py * cw + px);
          img.data[off] = lut[3 * u];
          img.data[off + 1] = lut[3 * u + 1];
          img.data[off + 2] = lut[3 * u + 2];
          img.data[off + 3] = 255;
        }
      }
      ctx.putImageData(img, 0, 0);

      const svg = d3.select(stage)
        .append("svg")
        .attr("viewBox", `0 0 ${W} ${H}`)
        .attr("preserveAspectRatio", "xMidYMid meet");

      const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);
      g.append("image")
        .attr("href", cvs.toDataURL())
        .attr("x", 0).attr("y", 0)
        .attr("width", iw).attr("height", ih)
        .attr("preserveAspectRatio", "none");

      const xScale = d3.scaleLinear().domain([xmin, xmax]).range([0, iw]);
      const yScale = d3.scaleLinear().domain([ymin, ymax]).range([ih, 0]);

      g.append("g")
        .attr("class", "cp-axis")
        .attr("transform", `translate(0,${ih})`)
        .call(d3.axisBottom(xScale).ticks(6));
      g.append("g")
        .attr("class", "cp-axis")
        .call(d3.axisLeft(yScale).ticks(6));

      if (panel.xlabel) {
        g.append("text")
          .attr("x", iw / 2).attr("y", ih + margin.bottom - 6)
          .attr("text-anchor", "middle")
          .attr("fill", "var(--cp-text)")
          .style("font-size", "11.5px")
          .text(panel.xlabel);
      }
      if (panel.ylabel) {
        g.append("text")
          .attr("transform", "rotate(-90)")
          .attr("x", -ih / 2).attr("y", -margin.left + 15)
          .attr("text-anchor", "middle")
          .attr("fill", "var(--cp-text)")
          .style("font-size", "11.5px")
          .text(panel.ylabel);
      }

      appendColorbar(svg, margin.left + iw + 14, margin.top, 12, ih, zmin, zmax, panel.colormap || "viridis", panel.cbarLabel || "");
      styleAxes(svg);

      const bisectXC = d3.bisector((d) => d).center;
      const bisectYC = d3.bisector((d) => d).center;
      g.append("rect")
        .attr("x", 0).attr("y", 0)
        .attr("width", iw).attr("height", ih)
        .attr("fill", "transparent")
        .style("cursor", "crosshair")
        .on("mousemove", function (evt) {
          const [mx, my] = d3.pointer(evt, this);
          const xv = xScale.invert(mx);
          const yv = yScale.invert(my);
          const ix = bisectXC(xArr, xv);
          const iy = bisectYC(yArr, yv);
          const val = zMat[iy][ix];
          showTooltip(
            evt,
            `<div class="canopy-d3-tooltip-head">${panel.title || "Heatmap"}</div>` +
            `<div class="canopy-d3-tooltip-row"><span>${panel.xlabel || "x"}</span><strong>${fmtNum(xArr[ix])}</strong></div>` +
            `<div class="canopy-d3-tooltip-row"><span>${panel.ylabel || "y"}</span><strong>${fmtNum(yArr[iy])}</strong></div>` +
            `<div class="canopy-d3-tooltip-row"><span>${panel.cbarLabel || "value"}</span><strong>${fmtNum(val)}</strong></div>`
          );
        })
        .on("mouseleave", hideTooltip);
    }

    draw();
    return draw;
  }

  function renderPolarPanel(panelEl, panel, rowHeight) {
    const stage = document.createElement("div");
    stage.className = "canopy-d3-stage";
    panelEl.appendChild(stage);

    function draw() {
      stage.innerHTML = "";
      const W = Math.max(280, stage.clientWidth || 420);
      const H = panel.height || rowHeight || 310;
      const cbSpace = 72;
      const pad = 26;
      const cx = (W - cbSpace) / 2;
      const cy = H / 2 + 2;
      const R = Math.max(55, Math.min(cx - pad, H / 2 - pad));

      const vzaArr = panel.vza; // [0, 5, ..., 80]
      const vazArr = panel.vaz; // [0, 5, ..., 355]
      const zMat = panel.z;     // [nvza][nvaz]
      const nvza = vzaArr.length;
      const nvaz = vazArr.length;
      const rMax = vzaArr[nvza - 1];
      const dvaz = vazArr.length > 1 ? vazArr[1] - vazArr[0] : 360;

      let zmin = Infinity, zmax = -Infinity;
      if (panel.clim) {
        zmin = panel.clim[0];
        zmax = panel.clim[1];
      } else {
        for (let ir = 0; ir < nvza; ir++) {
          for (let ia = 0; ia < nvaz; ia++) {
            const v = zMat[ir][ia];
            if (Number.isFinite(v)) {
              if (v < zmin) zmin = v;
              if (v > zmax) zmax = v;
            }
          }
        }
        if (!Number.isFinite(zmin)) { zmin = 0; zmax = 1; }
        if (zmin === zmax) { zmax = zmin + 1e-6; }
      }

      // Render smooth polar disc onto offscreen canvas
      const lut = buildColorLUT(panel.colormap || "viridis");
      const dim = Math.max(220, Math.round(2 * R * 1.35));
      const cvs = document.createElement("canvas");
      cvs.width = dim;
      cvs.height = dim;
      const ctx = cvs.getContext("2d");
      const img = ctx.createImageData(dim, dim);
      const bisectR = d3.bisector((d) => d).right;

      function samplePolar(vzaDeg, vazDeg) {
        let r1 = bisectR(vzaArr, vzaDeg);
        r1 = Math.max(1, Math.min(nvza - 1, r1));
        const r0 = r1 - 1;
        const tr = vzaArr[r1] === vzaArr[r0] ? 0 : Math.max(0, Math.min(1, (vzaDeg - vzaArr[r0]) / (vzaArr[r1] - vzaArr[r0])));
        const az = ((vazDeg % 360) + 360) % 360;
        const fAz = az / dvaz;
        const a0 = Math.floor(fAz) % nvaz;
        const a1 = (a0 + 1) % nvaz;
        const ta = fAz - Math.floor(fAz);
        const v00 = zMat[r0][a0], v01 = zMat[r0][a1];
        const v10 = zMat[r1][a0], v11 = zMat[r1][a1];
        return (1 - tr) * ((1 - ta) * v00 + ta * v01) + tr * ((1 - ta) * v10 + ta * v11);
      }

      const half = dim / 2;
      for (let py = 0; py < dim; py++) {
        const dy = py + 0.5 - half;
        for (let px = 0; px < dim; px++) {
          const dx = px + 0.5 - half;
          const rr = Math.hypot(dx, dy) / half;
          if (rr <= 1.0) {
            const vzaDeg = rr * rMax;
            // North = 0 deg (dx=0, dy<0), clockwise = +vaz (dx>0)
            const vazDeg = ((Math.atan2(dx, -dy) * 180) / Math.PI + 360) % 360;
            const val = samplePolar(vzaDeg, vazDeg);
            const u = Math.max(0, Math.min(255, Math.round(255 * (val - zmin) / (zmax - zmin))));
            const off = 4 * (py * dim + px);
            img.data[off] = lut[3 * u];
            img.data[off + 1] = lut[3 * u + 1];
            img.data[off + 2] = lut[3 * u + 2];
            img.data[off + 3] = 255;
          }
        }
      }
      ctx.putImageData(img, 0, 0);

      const svg = d3.select(stage)
        .append("svg")
        .attr("viewBox", `0 0 ${W} ${H}`)
        .attr("preserveAspectRatio", "xMidYMid meet");

      svg.append("image")
        .attr("href", cvs.toDataURL())
        .attr("x", cx - R).attr("y", cy - R)
        .attr("width", 2 * R).attr("height", 2 * R);

      const overlay = svg.append("g");
      // Concentric zenith circles
      const rTicks = [20, 40, 60, 80].filter((v) => v <= rMax);
      rTicks.forEach((rv) => {
        const rad = (rv / rMax) * R;
        overlay.append("circle")
          .attr("cx", cx).attr("cy", cy).attr("r", rad)
          .attr("fill", "none")
          .attr("stroke", rv === rMax ? "var(--cp-axis)" : "rgba(255,255,255,0.42)")
          .attr("stroke-width", rv === rMax ? 1.2 : 0.85)
          .attr("stroke-dasharray", rv === rMax ? null : "3,3");
        overlay.append("text")
          .attr("x", cx + 4)
          .attr("y", cy - rad + 10)
          .attr("fill", rv === rMax ? "var(--cp-text)" : "#ffffff")
          .attr("stroke", rv === rMax ? "none" : "rgba(0,0,0,0.55)")
          .attr("stroke-width", rv === rMax ? 0 : 2.2)
          .attr("paint-order", "stroke")
          .style("font-size", "9.5px")
          .text(rv + "°");
      });

      // Azimuth spokes & labels (0 deg at top, clockwise)
      for (let az = 0; az < 360; az += 45) {
        const ang = ((az - 90) * Math.PI) / 180;
        const x2 = cx + R * Math.cos(ang);
        const y2 = cy + R * Math.sin(ang);
        overlay.append("line")
          .attr("x1", cx).attr("y1", cy)
          .attr("x2", x2).attr("y2", y2)
          .attr("stroke", "rgba(255,255,255,0.35)")
          .attr("stroke-width", 0.8);
        const lx = cx + (R + 13) * Math.cos(ang);
        const ly = cy + (R + 13) * Math.sin(ang);
        overlay.append("text")
          .attr("x", lx)
          .attr("y", ly + 3.5)
          .attr("text-anchor", "middle")
          .attr("fill", "var(--cp-text)")
          .style("font-size", "10px")
          .text(az + "°");
      }

      // Hotspot / reference markers
      (panel.markers || []).forEach((mk) => {
        const rad = (mk.vza / rMax) * R;
        const ang = ((mk.vaz - 90) * Math.PI) / 180;
        const mx = cx + rad * Math.cos(ang);
        const my = cy + rad * Math.sin(ang);
        overlay.append("polygon")
          .attr("points", starPoints(mx, my, 7.2, 3.2))
          .attr("fill", mk.color || "#ef4444")
          .attr("stroke", "#ffffff")
          .attr("stroke-width", 1.2);
      });

      appendColorbar(
        svg,
        W - cbSpace + 12,
        cy - R + 8,
        12,
        Math.max(80, 2 * R - 16),
        zmin,
        zmax,
        panel.colormap || "viridis",
        panel.cbarLabel || ""
      );
      styleAxes(svg);

      svg.append("circle")
        .attr("cx", cx).attr("cy", cy).attr("r", R)
        .attr("fill", "transparent")
        .style("cursor", "crosshair")
        .on("mousemove", function (evt) {
          const [mx, my] = d3.pointer(evt, this);
          const dx = mx - cx;
          const dy = my - cy;
          const rr = Math.hypot(dx, dy) / R;
          if (rr > 1.02) {
            hideTooltip();
            return;
          }
          const vzaDeg = Math.min(rMax, rr * rMax);
          const vazDeg = ((Math.atan2(dx, -dy) * 180) / Math.PI + 360) % 360;
          const val = samplePolar(vzaDeg, vazDeg);
          showTooltip(
            evt,
            `<div class="canopy-d3-tooltip-head">${panel.title || "Polar BRDF"}</div>` +
            `<div class="canopy-d3-tooltip-row"><span>VZA</span><strong>${vzaDeg.toFixed(1)}°</strong></div>` +
            `<div class="canopy-d3-tooltip-row"><span>VAZ (RAA)</span><strong>${vazDeg.toFixed(1)}°</strong></div>` +
            `<div class="canopy-d3-tooltip-row"><span>${panel.cbarLabel || "Value"}</span><strong>${fmtNum(val)}</strong></div>`
          );
        })
        .on("mouseleave", hideTooltip);
    }

    draw();
    return draw;
  }

  function renderFigure(container, spec) {
    if (!window.d3) return;
    container.innerHTML = "";

    if (spec.title || spec.subtitle) {
      const hdr = document.createElement("div");
      hdr.className = "canopy-d3-header";
      if (spec.title) {
        const h = document.createElement("div");
        h.className = "canopy-d3-title";
        h.textContent = spec.title;
        hdr.appendChild(h);
      }
      if (spec.subtitle) {
        const sub = document.createElement("div");
        sub.className = "canopy-d3-subtitle";
        sub.textContent = spec.subtitle;
        hdr.appendChild(sub);
      }
      container.appendChild(hdr);
    }

    let currentFrame = 0;
    let playTimer = null;
    let frameLabelEl = null;
    let sliderEl = null;
    let playBtnEl = null;

    if (spec.frames && spec.frames.length > 0) {
      const ctrl = document.createElement("div");
      ctrl.className = "canopy-d3-controls";

      playBtnEl = document.createElement("button");
      playBtnEl.type = "button";
      playBtnEl.className = "canopy-d3-play-btn";
      playBtnEl.textContent = "▶ Play";

      sliderEl = document.createElement("input");
      sliderEl.type = "range";
      sliderEl.className = "canopy-d3-slider";
      sliderEl.min = "0";
      sliderEl.max = String(spec.frames.length - 1);
      sliderEl.step = "1";
      sliderEl.value = "0";

      frameLabelEl = document.createElement("span");
      frameLabelEl.className = "canopy-d3-frame-label";
      frameLabelEl.textContent = spec.frames[0].label || "";

      ctrl.appendChild(playBtnEl);
      ctrl.appendChild(sliderEl);
      ctrl.appendChild(frameLabelEl);
      container.appendChild(ctrl);
    }

    const grid = document.createElement("div");
    grid.className = "canopy-d3-grid";
    container.appendChild(grid);

    let panelDrawers = [];

    function renderPanels(panels) {
      grid.innerHTML = "";
      panelDrawers = [];
      const cols = spec.cols || Math.min(3, panels.length || 1);
      grid.style.setProperty("--cp-cols", String(cols));

      panels.forEach((p) => {
        const pEl = document.createElement("div");
        pEl.className = "canopy-d3-panel";
        if (p.title) {
          const pt = document.createElement("div");
          pt.className = "canopy-d3-panel-title";
          pt.textContent = p.title;
          pEl.appendChild(pt);
        }
        grid.appendChild(pEl);
        const kind = p.type || "cartesian";
        if (kind === "heatmap") {
          panelDrawers.push(renderHeatmapPanel(pEl, p, spec.rowHeight));
        } else if (kind === "polar") {
          panelDrawers.push(renderPolarPanel(pEl, p, spec.rowHeight));
        } else {
          panelDrawers.push(renderCartesianPanel(pEl, p, spec.rowHeight));
        }
      });
    }

    function applyFrame(idx) {
      currentFrame = idx;
      if (sliderEl) sliderEl.value = String(idx);
      const fr = spec.frames[idx];
      if (frameLabelEl) frameLabelEl.textContent = fr.label || "";
      renderPanels(fr.panels || []);
    }

    if (spec.frames && spec.frames.length > 0) {
      applyFrame(0);
      sliderEl.addEventListener("input", () => {
        if (playTimer) {
          clearInterval(playTimer);
          playTimer = null;
          playBtnEl.textContent = "▶ Play";
        }
        applyFrame(Number(sliderEl.value));
      });
      playBtnEl.addEventListener("click", () => {
        if (playTimer) {
          clearInterval(playTimer);
          playTimer = null;
          playBtnEl.textContent = "▶ Play";
        } else {
          playBtnEl.textContent = "⏸ Pause";
          playTimer = setInterval(() => {
            const next = (currentFrame + 1) % spec.frames.length;
            applyFrame(next);
          }, spec.intervalMs || 650);
        }
      });
    } else {
      renderPanels(spec.panels || []);
    }

    container.__canopyRedraw = function () {
      panelDrawers.forEach((fn) => fn());
    };
  }

  function resolveVendorD3Url() {
    const scripts = document.querySelectorAll("script[src]");
    for (let i = 0; i < scripts.length; i++) {
      const s = scripts[i].getAttribute("src") || "";
      if (s.indexOf("canopy_plots.js") !== -1) {
        return s.replace(/canopy_plots\.js(\?.*)?$/, "vendor/d3.v7.min.js");
      }
    }
    return "assets/vendor/d3.v7.min.js";
  }

  let d3Loading = false;
  const d3Callbacks = [];

  function ensureD3(cb) {
    if (window.d3 && window.d3.select) {
      cb();
      return;
    }
    d3Callbacks.push(cb);
    if (d3Loading) return;
    d3Loading = true;
    const d3Url = resolveVendorD3Url();

    function flush() {
      while (d3Callbacks.length > 0) {
        const fn = d3Callbacks.shift();
        fn();
      }
    }

    if (typeof window.require === "function" && window.define && window.define.amd) {
      window.require([d3Url], function (d3mod) {
        window.d3 = d3mod;
        flush();
      });
    } else {
      const tag = document.createElement("script");
      tag.src = d3Url;
      tag.onload = flush;
      document.head.appendChild(tag);
    }
  }

  function initAll() {
    const nodes = document.querySelectorAll(".canopy-d3-plot");
    if (nodes.length === 0) return;
    ensureD3(function () {
      nodes.forEach((el) => {
        if (el.dataset.cpInitialized === "1") return;
        const script = el.querySelector("script.canopy-plot-spec");
        if (!script) return;
        try {
          const spec = JSON.parse(script.textContent);
          el.dataset.cpInitialized = "1";
          renderFigure(el, spec);
        } catch (err) {
          console.error("CanopyPlots failed to parse figure spec:", err);
        }
      });
    });
  }

  window.CanopyPlots = { render: renderFigure, initAll: initAll };

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initAll);
  } else {
    initAll();
  }

  let resizeRaf = 0;
  window.addEventListener("resize", () => {
    if (resizeRaf) cancelAnimationFrame(resizeRaf);
    resizeRaf = requestAnimationFrame(() => {
      document.querySelectorAll(".canopy-d3-plot").forEach((el) => {
        if (typeof el.__canopyRedraw === "function") el.__canopyRedraw();
      });
    });
  });
})();
