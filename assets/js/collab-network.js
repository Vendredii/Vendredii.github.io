// Collaboration network on a map, in the style of VOSviewer.
//
//   circles  places from _data/collaborators.yml; institutions at the same
//            coordinates share one circle, sized by the number of papers
//            with an author there
//   lines    two places are linked when people from both are authors of the
//            same paper in _data/publications.yml (thicker = more papers)
//   colours  clusters of places that often publish together
//
// A co-author is placed by matching their name in publications.yml to the
// `people` lists in collaborators.yml. You are placed at the `home` entry.
(function () {
    var el = document.getElementById("collab-map");
    if (!el || !window.L) return;
    var data = window.CV_DATA;

    // Colours of VOSviewer's default cluster palette, in cluster-size order.
    var COLORS = ["#d62728", "#2ca02c", "#1f77b4", "#bcbd22", "#9467bd", "#17becf",
                  "#ff7f0e", "#8c564b", "#e377c2", "#7f7f7f", "#aec7e8", "#98df8a"];

    function clean(name) { return String(name).replace(/\*/g, "").trim(); }

    // ---- Places, grouped by location ----
    var groups = [], byKey = {}, placesOf = {};
    data.places.forEach(function (p) {
        if (typeof p.lat !== "number" || typeof p.lon !== "number") return;
        var key = p.lat.toFixed(2) + "," + p.lon.toFixed(2);
        if (!byKey[key]) {
            byKey[key] = { lat: p.lat, lon: p.lon, name: String(p.place || p.institution).split(",")[0],
                           place: p.place, home: false, items: [], papers: 0, links: {} };
            byKey[key].index = groups.length;
            groups.push(byKey[key]);
        }
        var g = byKey[key];
        g.items.push(p);
        if (p.home) g.home = true;
        String(p.people || "").split(",").map(clean).filter(Boolean).forEach(function (name) {
            placesOf[name] = placesOf[name] || [];
            if (placesOf[name].indexOf(g) < 0) placesOf[name].push(g);
        });
    });
    if (!groups.length) return;
    var home = groups.filter(function (g) { return g.home; })[0];
    if (home && data.me) placesOf[clean(data.me)] = [home];

    // ---- Links from shared papers ----
    (data.publications || []).forEach(function (pub) {
        var gs = [];
        (pub.authors || []).forEach(function (a) {
            (placesOf[clean(a)] || []).forEach(function (g) { if (gs.indexOf(g) < 0) gs.push(g); });
        });
        gs.forEach(function (g, i) {
            g.papers += 1;
            for (var j = i + 1; j < gs.length; j++) {
                g.links[gs[j].index] = (g.links[gs[j].index] || 0) + 1;
                gs[j].links[g.index] = (gs[j].links[g.index] || 0) + 1;
            }
        });
    });

    // ---- Clusters: modularity-based local moving, as VOSviewer does ----
    var strength = groups.map(function (g) {
        return Object.keys(g.links).reduce(function (s, j) { return s + g.links[j]; }, 0);
    });
    var total = strength.reduce(function (a, b) { return a + b; }, 0) || 1;
    var comm = groups.map(function (g, i) { return i; }), tot = strength.slice();
    for (var pass = 0, moved = true; moved && pass < 50; pass++) {
        moved = false;
        groups.forEach(function (g, i) {
            var old = comm[i], toComm = {};
            tot[old] -= strength[i];
            Object.keys(g.links).forEach(function (j) { toComm[comm[j]] = (toComm[comm[j]] || 0) + g.links[j]; });
            var best = old, bestGain = (toComm[old] || 0) - tot[old] * strength[i] / total;
            Object.keys(toComm).forEach(function (c) {
                var gain = toComm[c] - tot[c] * strength[i] / total;
                if (gain > bestGain + 1e-9) { best = +c; bestGain = gain; }
            });
            comm[i] = best; tot[best] += strength[i];
            if (best !== old) moved = true;
        });
    }
    var sizes = {};
    comm.forEach(function (c) { sizes[c] = (sizes[c] || 0) + 1; });
    var order = Object.keys(sizes).sort(function (a, b) { return sizes[b] - sizes[a] || a - b; });
    groups.forEach(function (g, i) { g.color = COLORS[order.indexOf(String(comm[i])) % COLORS.length]; });

    // ---- Map ----
    var map = L.map(el, { scrollWheelZoom: false, worldCopyJump: true, zoomSnap: 0.25 });
    // OpenStreetMap, turned grey in style.css so the coloured network stands out.
    L.tileLayer("https://tile.openstreetmap.org/{z}/{x}/{y}.png", {
        maxZoom: 18,
        attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
    }).addTo(map);

    function radius(g) { return 5 + 3 * Math.sqrt(g.papers); }
    function fontSize(g) { return 11 + 1.6 * Math.sqrt(g.papers); }

    // Curved links, bending the same way as in VOSviewer.
    var lines = [];
    groups.forEach(function (a) {
        Object.keys(a.links).forEach(function (j) {
            var b = groups[j];
            if (b.index < a.index) return;
            var p = a.lon <= b.lon ? a : b, q = p === a ? b : a;
            var cx = (p.lon + q.lon) / 2 - (q.lat - p.lat) * 0.2, cy = (p.lat + q.lat) / 2 + (q.lon - p.lon) * 0.2;
            var pts = [];
            for (var t = 0; t <= 1.0001; t += 1 / 32) {
                pts.push([(1 - t) * (1 - t) * p.lat + 2 * (1 - t) * t * cy + t * t * q.lat,
                          (1 - t) * (1 - t) * p.lon + 2 * (1 - t) * t * cx + t * t * q.lon]);
            }
            var w = a.links[j];
            var line = L.polyline(pts, {
                color: (a.papers >= b.papers ? a : b).color, weight: 0.8 + 1.2 * w, opacity: 0.45, interactive: false
            }).addTo(map);
            line.ends = [a, b];
            lines.push(line);
        });
    });

    var focus = null;
    groups.forEach(function (g) {
        var popup = document.createElement("div");
        var head = document.createElement("strong");
        head.textContent = g.place + " · " + g.papers + (g.papers === 1 ? " paper" : " papers");
        popup.appendChild(head);
        g.items.forEach(function (p) {
            popup.appendChild(document.createElement("hr"));
            var name = document.createElement("div");
            name.className = "popup-inst";
            name.textContent = p.institution;
            popup.appendChild(name);
            if (p.people) popup.appendChild(document.createTextNode(p.people));
        });

        g.marker = L.circleMarker([g.lat, g.lon], {
            radius: radius(g), color: "#fff", weight: 1.5, fillColor: g.color, fillOpacity: 0.9
        }).bindPopup(popup, { maxHeight: 260 }).addTo(map);
        g.marker.on("mouseover", function () { highlight(g); });
        g.marker.on("mouseout", function () { highlight(null); });

        g.label = L.tooltip({ permanent: true, direction: "center", className: "net-label", interactive: false })
            .setLatLng([g.lat, g.lon]).setContent(g.name).addTo(map);
    });

    // Show as many labels as fit without overlapping, biggest places first.
    var byPapers = groups.slice().sort(function (a, b) { return b.papers - a.papers || a.name.localeCompare(b.name); });
    function placeLabels() {
        var placed = [];
        byPapers.forEach(function (g) {
            var pt = map.latLngToContainerPoint([g.lat, g.lon]);
            var fs = fontSize(g), w = g.name.length * fs * 0.6, h = fs * 1.2;
            var box = [pt.x - w / 2, pt.y - h / 2, pt.x + w / 2, pt.y + h / 2];
            var near = focus && (g === focus || focus.links[g.index]);
            var clash = placed.some(function (b) { return box[0] < b[2] && box[2] > b[0] && box[1] < b[3] && box[3] > b[1]; });
            var show = focus ? near : !clash;
            if (show) placed.push(box);
            var node = g.label.getElement();
            if (node) {
                node.style.display = show ? "" : "none";
                node.style.fontSize = fs + "px";
            }
        });
    }

    // Hovering a place highlights it, its partners and the links between them.
    function highlight(g) {
        focus = g;
        groups.forEach(function (h) {
            var on = !g || h === g || g.links[h.index];
            h.marker.setStyle({ fillOpacity: on ? 0.9 : 0.15, opacity: on ? 1 : 0.15 });
        });
        lines.forEach(function (line) {
            var on = !g || line.ends.indexOf(g) >= 0;
            line.setStyle({ opacity: g ? (on ? 0.85 : 0.05) : 0.45 });
        });
        placeLabels();
    }

    map.fitBounds(groups.map(function (g) { return [g.lat, g.lon]; }), { padding: [40, 40], maxZoom: 5 });
    map.on("zoomend moveend", placeLabels);
    placeLabels();
    // Draw the biggest places on top.
    byPapers.slice().reverse().forEach(function (g) { g.marker.bringToFront(); });
})();
