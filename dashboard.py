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

GITHUB_SVG = (
    '<svg xmlns="http://www.w3.org/2000/svg" width="14" height="14" viewBox="0 0 24 24" '
    'fill="currentColor" style="vertical-align:-2px;">'
    '<path d="M12 0C5.37 0 0 5.37 0 12c0 5.3 3.438 9.8 8.205 11.385.6.113.82-.258.82-.577'
    " 0-.285-.01-1.04-.015-2.04-3.338.724-4.042-1.61-4.042-1.61C4.422 18.07 3.633 17.7"
    " 3.633 17.7c-1.087-.744.084-.729.084-.729 1.205.084 1.838 1.236 1.838 1.236"
    " 1.07 1.835 2.809 1.305 3.495.998.108-.776.417-1.305.76-1.605-2.665-.3-5.466-1.332"
    "-5.466-5.93 0-1.31.465-2.38 1.235-3.22-.135-.303-.54-1.523.105-3.176 0 0"
    " 1.005-.322 3.3 1.23.96-.267 1.98-.399 3-.405 1.02.006 2.04.138 3 .405"
    " 2.28-1.552 3.285-1.23 3.285-1.23.645 1.653.24 2.873.12 3.176.765.84"
    " 1.23 1.91 1.23 3.22 0 4.61-2.805 5.625-5.475 5.92.42.36.81 1.096.81"
    " 2.22 0 1.606-.015 2.896-.015 3.286 0 .315.21.69.825.57C20.565 21.795"
    ' 24 17.295 24 12c0-6.63-5.37-12-12-12"/></svg>'
)

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

.tool-name a {
    color: inherit;
    text-decoration: none;
    border-bottom: 1px dotted transparent;
    transition: color 0.15s, border-color 0.15s;
}
.tool-name a:hover {
    color: var(--ac);
    border-bottom-color: var(--ac);
}

.card-links { display: flex; gap: 8px; flex-wrap: wrap; margin-top: 0.55rem; }
.card-link {
    display: inline-flex; align-items: center; gap: 5px;
    font-size: 0.78rem; font-weight: 500;
    padding: 3px 10px; border-radius: 6px; text-decoration: none;
    border: 1px solid #e2e8f0; background: #f8fafc; color: #334155;
}
.card-link:hover { background: #f1f5f9; border-color: #cbd5e1; }
.card-link.gh    { background: #f6f8fa; color: #24292f; border-color: #d0d7de; }

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


def _extract_doi(text: str) -> str:
    """Return the first DOI URL found, or empty string."""
    m = re.search(r'https?://(?:dx\.)?doi\.org/10\.\d{4,}/[^\s\)\"\]]+', text)
    return m.group(0).rstrip('.,;') if m else ""


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
                'doi':        _extract_doi(block),
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

    # Title: link to paper DOI if available, otherwise plain text
    if tool['doi']:
        name_html = (f'<a href="{html.escape(tool["doi"])}" target="_blank" '
                     f'title="Open paper">{name}</a>')
    else:
        name_html = name

    # Inline buttons inside the card
    gh_btn = ""
    if tool['github_url']:
        gh_btn = (f'<a class="card-link gh" href="{html.escape(tool["github_url"])}" '
                  f'target="_blank">{GITHUB_SVG}&nbsp;GitHub</a>')

    proj_btn = ""
    if tool['url'] and tool['url'] != tool['github_url']:
        proj_btn = (f'<a class="card-link" href="{html.escape(tool["url"])}" '
                    f'target="_blank">🔗 Project</a>')

    st.markdown(f"""
<div class="tool-card" style="--ac:{color};">
  <div class="tool-name">{name_html}</div>
  <div class="chips">{cat_chip}{lang_chips}{year_chip}</div>
  <div class="tool-desc">{desc}</div>
  <div class="card-links">{gh_btn}{proj_btn}</div>
</div>""", unsafe_allow_html=True)


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
