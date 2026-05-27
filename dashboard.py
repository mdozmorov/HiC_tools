#!/usr/bin/env python3
"""HiCatalog — Interactive Hi-C Tools Explorer"""

import html
import re
import streamlit as st
import pandas as pd
from collections import defaultdict
from pathlib import Path

MAX_DESCRIPTION_LENGTH = 280

CATEGORY_COLORS = [
    "#2563eb", "#16a34a", "#9333ea", "#dc2626", "#0891b2",
    "#d97706", "#db2777", "#059669", "#7c3aed", "#ea580c",
    "#0284c7", "#65a30d", "#c026d3", "#e11d48", "#0d9488",
]

LANGUAGES = [
    "Python", "R", "Java", "C++", "Perl", "Julia",
    "Matlab", "JavaScript", "Nextflow", "Snakemake", "Shell",
]

st.set_page_config(
    page_title="HiCatalog",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
html, body, [class*="css"] { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }

.hero {
    background: linear-gradient(135deg, #0f172a 0%, #1e3a5f 60%, #0c4a6e 100%);
    color: white;
    padding: 2rem 2.5rem;
    border-radius: 14px;
    margin-bottom: 1.75rem;
}
.hero h1 { color: white; font-size: 1.9rem; font-weight: 700; margin: 0 0 0.4rem 0; }
.hero p  { color: rgba(255,255,255,0.7); margin: 0; font-size: 0.9rem; }

.stats-row { display: flex; gap: 1rem; margin-bottom: 1.5rem; }
.stat-card {
    flex: 1;
    background: white;
    border: 1px solid #e2e8f0;
    border-radius: 10px;
    padding: 1rem 1.25rem;
    text-align: center;
    box-shadow: 0 1px 4px rgba(0,0,0,0.05);
}
.stat-card .num   { font-size: 1.85rem; font-weight: 700; color: #0f172a; }
.stat-card .lbl   { font-size: 0.72rem; color: #94a3b8; text-transform: uppercase;
                    letter-spacing: 0.06em; margin-top: 2px; }

.tool-card {
    background: white;
    border: 1px solid #e2e8f0;
    border-left: 4px solid var(--ac, #2563eb);
    border-radius: 10px;
    padding: 1.1rem 1.4rem;
    margin: 0.6rem 0;
    box-shadow: 0 1px 4px rgba(0,0,0,0.04);
}
.tool-card:hover { box-shadow: 0 4px 16px rgba(0,0,0,0.09); }

.tool-name { font-size: 1rem; font-weight: 600; color: #0f172a; margin-bottom: 0.3rem; }
.tool-desc { font-size: 0.85rem; color: #475569; line-height: 1.65; margin: 0.4rem 0 0.6rem; }

.chips { margin-bottom: 0.5rem; }
.chip {
    display: inline-block;
    padding: 2px 9px;
    border-radius: 99px;
    font-size: 0.7rem;
    font-weight: 500;
    margin: 2px 3px 2px 0;
}
.chip-cat  { background: var(--ac-bg, #eff6ff); color: var(--ac, #2563eb); }
.chip-lang { background: #f0fdf4; color: #15803d; }
.chip-year { background: #fff7ed; color: #c2410c; }

/* Pull the native-widget button row up so it sits flush under the card */
.tool-card { padding-bottom: 0.5rem; margin-bottom: 0; }
.tool-card + div[data-testid="stHorizontalBlock"] {
    margin-top: -0.4rem;
    margin-bottom: 0.6rem;
    padding: 0 0 0 1.4rem;
}

div[data-testid="stSidebar"] { background: #f8fafc; }
.block-container { padding-top: 1.25rem; }
section[data-testid="stSidebar"] > div { padding-top: 1rem; }
</style>
""", unsafe_allow_html=True)


def _extract_languages(text: str) -> list[str]:
    tl = text.lower()
    return [l for l in LANGUAGES if l.lower() in tl][:3]


def _extract_year(text: str) -> int | None:
    hits = re.findall(r'\b(20[0-9]{2}|199[0-9])\b', text)
    return max(int(y) for y in hits) if hits else None


@st.cache_data
def parse_readme(_mtime: float) -> tuple[list, dict]:
    """Parse README.md. `_mtime` is the file's modification time — passing it
    as an argument makes the cache invalidate when README.md is edited."""
    readme_path = Path(__file__).parent / "README.md"
    try:
        content = readme_path.read_text(encoding="utf-8")
    except OSError as e:
        st.error(f"Could not read README.md: {e}")
        return [], {}

    tools: list[dict] = []
    categories: dict = defaultdict(list)
    sections = re.split(r'\n## ', content)
    current_category = "General"

    SKIP = {"Table of content", "How to View This README", "Interactive Dashboard"}

    for section in sections:
        lines = section.split('\n')
        if not lines:
            continue
        m = re.match(r'^([^#\n]+)', lines[0])
        if m:
            current_category = m.group(1).strip()
        if current_category in SKIP:
            continue

        i = 0
        while i < len(lines):
            tm = re.match(r'^-\s*(?:<a name="([^"]+)">)?\[([^\]]+)\]\(([^)]+)\)',
                          lines[i])
            if not tm:
                i += 1
                continue

            # Collect the full block for this entry
            j = i + 1
            while j < len(lines) and not re.match(r'^-\s*(?:<a name=|)\[', lines[j]):
                j += 1
            block = "\n".join(lines[i:j])

            tool_id  = tm.group(1) or ""
            tool_name = tm.group(2)
            tool_url  = tm.group(3)

            # Description: text on the opening line, before <details>
            desc = lines[i][tm.end():].strip()
            desc = re.sub(r'<details>.*', '', desc, flags=re.DOTALL).strip()
            if len(desc) > MAX_DESCRIPTION_LENGTH:
                desc = desc[:MAX_DESCRIPTION_LENGTH] + "…"

            # Citation: full contents of <details>…</details>, minus <summary>
            citation = ""
            dm = re.search(r'<details>(.*?)</details>', block, re.DOTALL)
            if dm:
                inner = re.sub(r'<summary>.*?</summary>', '', dm.group(1),
                               flags=re.DOTALL)
                inner = re.sub(r'<[^>]+>', '', inner)
                citation = ' '.join(inner.split()).strip()
                if len(citation) > 600:
                    citation = citation[:600] + "…"

            # GitHub URL
            github_url = ""
            if "github.com" in tool_url:
                github_url = tool_url
            else:
                gh = re.search(r'https://github\.com/[^\s\)\"\]]+', block)
                if gh:
                    github_url = gh.group(0).rstrip('.,;')

            entry = {
                'name':       tool_name,
                'url':        tool_url,
                'github_url': github_url,
                'category':   current_category,
                'description': desc,
                'citation':   citation,
                'year':       _extract_year(block),
                'languages':  _extract_languages(block),
                'id':         tool_id,
            }
            tools.append(entry)
            categories[current_category].append(entry)
            i = j

    return tools, dict(categories)


def _color(category: str, all_cats: list[str]) -> str:
    idx = sorted(all_cats).index(category) % len(CATEGORY_COLORS)
    return CATEGORY_COLORS[idx]


def _card(tool: dict, color: str) -> None:
    name  = html.escape(tool['name'])
    desc  = html.escape(tool['description'])
    ac_bg = color + "18"   # ~10% opacity tint

    # chips
    cat_chip  = (f'<span class="chip chip-cat" style="--ac:{color};--ac-bg:{ac_bg};">'
                 f'{html.escape(tool["category"])}</span>')
    lang_chips = "".join(f'<span class="chip chip-lang">{l}</span>'
                         for l in tool['languages'])
    year_chip  = (f'<span class="chip chip-year">{tool["year"]}</span>'
                  if tool['year'] else "")

    st.markdown(f"""
<div class="tool-card" style="--ac:{color};">
  <div class="tool-name">{name}</div>
  <div class="chips">{cat_chip}{lang_chips}{year_chip}</div>
  <div class="tool-desc">{desc}</div>
</div>""", unsafe_allow_html=True)

    has_gh   = bool(tool['github_url'])
    has_proj = bool(tool['url']) and tool['url'] != tool['github_url']
    has_cite = bool(tool['citation'])

    if has_gh or has_proj or has_cite:
        cols = st.columns([1, 1, 1, 6])
        idx = 0
        if has_gh:
            cols[idx].link_button("🐙 GitHub", tool['github_url'],
                                  use_container_width=True)
            idx += 1
        if has_proj:
            cols[idx].link_button("🔗 Project", tool['url'],
                                  use_container_width=True)
            idx += 1
        if has_cite:
            with cols[idx].popover("📋 Citation", use_container_width=True):
                st.code(tool['citation'], language=None)


def main():
    st.markdown("""
<div class="hero">
  <h1>🧬 HiCatalog</h1>
  <p>The searchable encyclopedia of Hi-C data analysis tools</p>
</div>""", unsafe_allow_html=True)

    readme_mtime = (Path(__file__).parent / "README.md").stat().st_mtime
    tools, categories = parse_readme(readme_mtime)
    if not tools:
        st.error("README.md not found or no tools parsed.")
        return

    all_cats   = list(categories.keys())
    color_map  = {c: _color(c, all_cats) for c in all_cats}

    # ── Sidebar ──────────────────────────────────────────────────────────────
    st.sidebar.markdown("## Filter")
    search   = st.sidebar.text_input("Search", placeholder="Name, keyword, language…")
    sel_cat  = st.sidebar.selectbox("Category", ["All"] + sorted(categories.keys()))

    years = sorted({t['year'] for t in tools if t['year']})
    if years:
        yr = st.sidebar.slider("Publication year",
                                min_value=years[0], max_value=years[-1],
                                value=(years[0], years[-1]))
    else:
        yr = None

    all_langs   = sorted({l for t in tools for l in t['languages']})
    sel_langs   = st.sidebar.multiselect("Language", all_langs)

    sort_by   = st.sidebar.selectbox(
        "Sort by",
        ["Newest first", "Oldest first", "A → Z"],
    )
    view_mode   = st.sidebar.radio("View", ["Cards", "Table"])

    st.sidebar.markdown("---")
    st.sidebar.caption(f"{len(tools)} tools · {len(categories)} categories")

    # ── Filter ────────────────────────────────────────────────────────────────
    q = search.lower()
    filtered = [
        t for t in tools
        if (not q or q in t['name'].lower() or q in t['description'].lower()
                   or q in t['category'].lower()
                   or any(q in l.lower() for l in t['languages']))
        and (sel_cat == "All" or t['category'] == sel_cat)
        and (not yr or not t['year'] or yr[0] <= t['year'] <= yr[1])
        and (not sel_langs or any(l in t['languages'] for l in sel_langs))
    ]

    if sort_by == "Newest first":
        filtered.sort(key=lambda t: t['year'] or 0, reverse=True)
    elif sort_by == "Oldest first":
        filtered.sort(key=lambda t: t['year'] or 9999)
    else:
        filtered.sort(key=lambda t: t['name'].lower())

    # ── Stats ─────────────────────────────────────────────────────────────────
    c1, c2, c3 = st.columns(3)
    for col, num, lbl in zip(
        [c1, c2, c3],
        [len(tools), len(categories), len(filtered)],
        ["Total Tools", "Categories", "Showing"],
    ):
        with col:
            st.markdown(
                f'<div class="stat-card"><div class="num">{num}</div>'
                f'<div class="lbl">{lbl}</div></div>',
                unsafe_allow_html=True,
            )

    st.markdown("<br>", unsafe_allow_html=True)

    if not filtered:
        st.info("No tools match the current filters.")
        return

    # ── Results ───────────────────────────────────────────────────────────────
    if view_mode == "Cards":
        for t in filtered:
            _card(t, color_map.get(t['category'], CATEGORY_COLORS[0]))

    else:
        df = pd.DataFrame([{
            "Name":        t['name'],
            "Category":    t['category'],
            "Year":        t['year'] or "",
            "Languages":   ", ".join(t['languages']),
            "Description": t['description'],
            "Link":        t['url'],
            "GitHub":      t['github_url'],
        } for t in filtered])
        st.dataframe(
            df,
            column_config={
                "Link":   st.column_config.LinkColumn("Link"),
                "GitHub": st.column_config.LinkColumn("GitHub"),
            },
            hide_index=True,
            width="stretch",
        )

    st.markdown("""
<div style="text-align:center;color:#94a3b8;padding:1.5rem 0 0.5rem;font-size:0.78rem;">
  <a href="https://github.com/mdozmorov/HiC_tools" target="_blank"
     style="color:#94a3b8;">mdozmorov/HiC_tools</a>
  &nbsp;·&nbsp;
  <a href="https://github.com/mdozmorov/HiC_tools/blob/master/CONTRIBUTING.md"
     target="_blank" style="color:#94a3b8;">Contribute</a>
</div>""", unsafe_allow_html=True)


if __name__ == "__main__":
    main()
