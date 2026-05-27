#!/usr/bin/env python3
"""
Hi-C Tools Interactive Dashboard
A Streamlit-based interactive visualization for exploring Hi-C data analysis tools
"""

import streamlit as st
import pandas as pd
import re
from collections import defaultdict
from pathlib import Path

MAX_DESCRIPTION_LENGTH = 300
MAX_KANBAN_CATEGORIES = 6
MAX_TOOLS_PER_COLUMN = 10

st.set_page_config(
    page_title="Hi-C Tools Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.markdown("""
<style>
    .tool-card {
        background-color: #f0f2f6;
        border-radius: 10px;
        padding: 20px;
        margin: 10px 0;
        border-left: 5px solid #4CAF50;
    }
    .category-badge {
        background-color: #4CAF50;
        color: white;
        padding: 5px 10px;
        border-radius: 15px;
        font-size: 12px;
        margin-right: 5px;
        display: inline-block;
    }
    .stat-box {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 20px;
        border-radius: 10px;
        text-align: center;
        margin: 10px;
    }
    .stat-number { font-size: 36px; font-weight: bold; }
    .stat-label { font-size: 14px; opacity: 0.9; }
    h1 { color: #2E86AB; }
    h2 { color: #5C8AB8; }
</style>
""", unsafe_allow_html=True)


@st.cache_data
def parse_readme():
    """Parse README.md and extract tool entries by section."""
    tools = []
    categories = defaultdict(list)

    readme_path = Path(__file__).parent / "README.md"
    try:
        content = readme_path.read_text(encoding="utf-8")
    except OSError as e:
        st.error(f"Could not read README.md: {e}")
        return [], {}

    sections = re.split(r'\n## ', content)
    current_category = "General"

    for section in sections:
        lines = section.split('\n')
        if not lines:
            continue

        category_match = re.match(r'^([^#\n]+)', lines[0])
        if category_match:
            current_category = category_match.group(1).strip()

        skip = {"Table of content", "How to View This README",
                "Interactive Dashboard"}
        if current_category in skip:
            continue

        for line in lines:
            tool_match = re.match(
                r'^-\s*(?:<a name="([^"]+)">)?\[([^\]]+)\]\(([^)]+)\)', line
            )
            if not tool_match:
                continue

            tool_id = tool_match.group(1) or ""
            tool_name = tool_match.group(2)
            tool_url = tool_match.group(3)

            description = line[tool_match.end():].strip()
            description = re.sub(r'<details>.*', '', description, flags=re.DOTALL)
            if len(description) > MAX_DESCRIPTION_LENGTH:
                description = description[:MAX_DESCRIPTION_LENGTH] + "..."

            tool_info = {
                'name': tool_name,
                'url': tool_url,
                'category': current_category,
                'description': description,
                'id': tool_id,
            }
            tools.append(tool_info)
            categories[current_category].append(tool_info)

    return tools, dict(categories)


def display_tool_card(tool):
    st.markdown(f"""
    <div class="tool-card">
        <h3>🔧 {tool['name']}</h3>
        <span class="category-badge">{tool['category']}</span>
        <p>{tool['description']}</p>
        <a href="{tool['url']}" target="_blank">🔗 Visit Project</a>
    </div>
    """, unsafe_allow_html=True)


def main():
    st.markdown("# 🧬 Hi-C Tools Explorer")

    img_path = Path(__file__).parent / "img" / "3C_technologies.png"
    if img_path.exists():
        st.image(str(img_path), caption="3C Technologies Overview",
                 use_container_width=True)

    st.markdown("---")

    tools, categories = parse_readme()

    if not tools:
        st.error("No tools found. Please check that README.md is present.")
        return

    st.sidebar.title("🔍 Filter & Search")
    search_query = st.sidebar.text_input("Search tools", "",
                                         placeholder="Type tool name or keyword...")
    all_categories = ["All"] + sorted(categories.keys())
    selected_category = st.sidebar.selectbox("Filter by Category", all_categories)
    view_mode = st.sidebar.radio("View Mode", ["Cards", "Table", "Kanban Board"])

    st.sidebar.markdown("---")
    st.sidebar.markdown("### 📊 Statistics")
    st.sidebar.metric("Total Tools", len(tools))
    st.sidebar.metric("Categories", len(categories))

    col1, col2, col3 = st.columns(3)
    filtered_count = (len(categories.get(selected_category, []))
                      if selected_category != "All" else len(tools))
    for col, number, label in zip(
        [col1, col2, col3],
        [len(tools), len(categories), filtered_count],
        ["Total Tools", "Categories", "Showing"],
    ):
        with col:
            st.markdown(f"""
            <div class="stat-box">
                <div class="stat-number">{number}</div>
                <div class="stat-label">{label}</div>
            </div>
            """, unsafe_allow_html=True)

    st.markdown("---")

    filtered_tools = [
        t for t in tools
        if (not search_query or
            search_query.lower() in t['name'].lower() or
            search_query.lower() in t['description'].lower() or
            search_query.lower() in t['category'].lower())
        and (selected_category == "All" or t['category'] == selected_category)
    ]

    if not filtered_tools:
        st.warning("No tools match your search criteria.")
        return

    st.markdown(f"## Showing {len(filtered_tools)} tool(s)")

    if view_mode == "Cards":
        for tool in filtered_tools:
            display_tool_card(tool)

    elif view_mode == "Table":
        df = pd.DataFrame(filtered_tools)[['name', 'category', 'description', 'url']]
        st.dataframe(
            df,
            column_config={
                "url": st.column_config.LinkColumn("Link"),
                "name": "Tool Name",
                "category": "Category",
                "description": "Description",
            },
            hide_index=True,
            use_container_width=True,
        )

    elif view_mode == "Kanban Board":
        kanban_cats = (
            list(dict.fromkeys(t['category'] for t in filtered_tools))[:MAX_KANBAN_CATEGORIES]
            if selected_category == "All"
            else [selected_category]
        )
        cols = st.columns(len(kanban_cats))
        for col, category in zip(cols, kanban_cats):
            with col:
                st.markdown(f"### {category}")
                for tool in [t for t in filtered_tools
                             if t['category'] == category][:MAX_TOOLS_PER_COLUMN]:
                    desc = tool['description'][:100]
                    if len(tool['description']) > 100:
                        desc += "..."
                    st.markdown(f"""
                    <div style="background:#e3f2fd;padding:10px;margin:5px 0;
                                border-radius:5px;border-left:3px solid #2196F3;">
                        <strong>{tool['name']}</strong><br/>
                        <small>{desc}</small><br/>
                        <a href="{tool['url']}" target="_blank">🔗</a>
                    </div>
                    """, unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("""
    <div style="text-align:center;color:#666;padding:20px;">
        <p>💡 Use the sidebar to filter tools by category or search for specific functionality.</p>
        <p>📚 Full list in <a href="https://github.com/mdozmorov/HiC_tools/blob/master/README.md" target="_blank">README.md</a>
        · 🌟 <a href="https://github.com/mdozmorov/HiC_tools/blob/master/CONTRIBUTING.md" target="_blank">Contribute</a></p>
    </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
