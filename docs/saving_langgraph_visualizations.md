# Saving LangGraph Workflow Visualizations

**Author:** Samuel Ahuno (ekwame001@gmail.com)
**Date:** 2025-12-24
**Purpose:** Reference guide for generating and saving LangGraph workflow diagrams

---

## Overview

LangGraph workflows can be visualized as flowcharts showing nodes, edges, and routing logic. This guide covers three methods to save these visualizations for documentation and inspection.

## Prerequisites

```bash
# Required packages
pip install langgraph langchain-core
```

## Method 1: Command Line Interface (Fastest)

### Usage

```bash
# Basic usage - saves to current directory
python your_agent.py --save-graph output_name

# Save to specific directory
python your_agent.py --save-graph ./docs/workflow_diagram

# Save to parent directory
python your_agent.py --save-graph ../graphs/my_workflow
```

### Output Files

- `output_name.png` - Visual diagram (typically 50-100KB)
- `output_name.mermaid` - Mermaid source code (typically 1-3KB)

### Example

```bash
# For the HITL workflow
/Users/ahunos/anaconda3/envs/pagents/bin/python blockA_single_Agent_HITL.py --save-graph hitl_workflow

# Output:
# Graph PNG saved to: hitl_workflow.png
# Mermaid code saved to: hitl_workflow.mermaid
```

---

## Method 2: Python Script or Notebook

### Implementation

Add this helper function to your agent script:

```python
def save_graph_visualization(graph, output_path: str = "workflow_graph"):
    """
    Save LangGraph visualization to files.

    Args:
        graph: Compiled LangGraph
        output_path: Base path for output files (without extension)

    Saves:
        - {output_path}.png: Graph visualization as PNG
        - {output_path}.mermaid: Mermaid diagram code as text
    """
    import os

    # Ensure output directory exists
    output_dir = os.path.dirname(output_path) or "."
    os.makedirs(output_dir, exist_ok=True)

    # Save PNG image
    png_data = graph.get_graph(xray=True).draw_mermaid_png()
    png_path = f"{output_path}.png"
    with open(png_path, "wb") as f:
        f.write(png_data)
    print(f"Graph PNG saved to: {png_path}")

    # Save Mermaid code as text
    mermaid_code = graph.get_graph(xray=True).draw_mermaid()
    mermaid_path = f"{output_path}.mermaid"
    with open(mermaid_path, "w") as f:
        f.write(mermaid_code)
    print(f"Mermaid code saved to: {mermaid_path}")

    return png_path, mermaid_path
```

### Usage in Python

```python
from your_agent import create_workflow_graph, save_graph_visualization

# Create and compile your graph
builder = create_workflow_graph()
graph = builder.compile()

# Save visualization
save_graph_visualization(graph, "my_workflow_diagram")
```

### Example

```python
from blockA_single_Agent_HITL import create_hitl_graph, save_graph_visualization

builder = create_hitl_graph()
graph = builder.compile()
save_graph_visualization(graph, "./docs/hitl_workflow")
```

---

## Method 3: Jupyter Notebook (Interactive Display)

### Display Inline

```python
from your_agent import create_workflow_graph
from IPython.display import Image, display

# Create and compile graph
builder = create_workflow_graph()
graph = builder.compile()

# Display inline in notebook
display(Image(graph.get_graph(xray=True).draw_mermaid_png()))
```

### Save from Notebook

```python
# Option A: Use the helper function
from your_agent import save_graph_visualization
save_graph_visualization(graph, "notebook_workflow")

# Option B: Manual save
png_data = graph.get_graph(xray=True).draw_mermaid_png()
with open("workflow.png", "wb") as f:
    f.write(png_data)

mermaid_code = graph.get_graph(xray=True).draw_mermaid()
with open("workflow.mermaid", "w") as f:
    f.write(mermaid_code)
```

### Example Notebook Cell

```python
# Cell 1: Import and create graph
from blockA_single_Agent_HITL import create_hitl_graph
from IPython.display import Image, display

builder = create_hitl_graph()
graph = builder.compile()

# Cell 2: Display
display(Image(graph.get_graph(xray=True).draw_mermaid_png()))

# Cell 3: Save to file
from blockA_single_Agent_HITL import save_graph_visualization
save_graph_visualization(graph, "./outputs/hitl_diagram")
```

---

## Working with Output Files

### PNG Images

**View:**
- macOS: `open workflow.png`
- Linux: `xdg-open workflow.png`
- Windows: `start workflow.png`
- Any OS: Double-click the file

**Use in documentation:**
```markdown
![Workflow Diagram](./workflow.png)
```

### Mermaid Files

**View online:**
1. Copy contents of `.mermaid` file
2. Visit https://mermaid.live
3. Paste code in editor

**View in GitHub:**
- GitHub automatically renders `.mermaid` files
- Create `workflow.md` with:
  ```markdown
  # Workflow Diagram

  ```mermaid
  [paste mermaid code here]
  ```
  ```

**View in VS Code:**
- Install "Mermaid Preview" extension
- Open `.mermaid` file
- Press `Cmd+Shift+P` (macOS) or `Ctrl+Shift+P` (Windows/Linux)
- Type "Mermaid: Preview"

---

## Adding CLI Support to Your Agent

### Step 1: Add import

```python
import argparse
```

### Step 2: Add argument parser

```python
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Your Agent Description")

    parser.add_argument(
        "--save-graph",
        type=str,
        metavar="PATH",
        help="Save graph visualization to specified path (without extension)"
    )

    args = parser.parse_args()

    # Handle graph saving
    if args.save_graph:
        builder = create_your_graph()
        graph = builder.compile()
        save_graph_visualization(graph, args.save_graph)
        print("Graph visualization saved. Exiting.")
        exit(0)

    # ... rest of your main logic
```

### Step 3: Test

```bash
python your_agent.py --save-graph test_output
```

---

## Understanding the Graph Visualization

### Graph Elements

**Nodes (Rectangles):**
- Represent functions/operations in your workflow
- Example: `detect_ambiguity`, `search_pubmed`, `validate_results`

**Edges (Lines):**
- **Solid lines** (`→`): Fixed edges (always follow this path)
- **Dotted lines** (`-.->` or `- . label . ->`): Conditional edges (routing based on state)

**Special Nodes:**
- `__start__`: Entry point
- `__end__`: Exit point(s)

### xray=True Parameter

```python
graph.get_graph(xray=True).draw_mermaid_png()
```

**Purpose:** Shows internal structure including:
- Subgraph internals
- All conditional routing options
- Hidden nodes

**When to use:**
- ✅ **Use `xray=True`** for debugging and understanding full workflow
- ❌ **Use `xray=False`** for simplified diagrams in documentation

---

## Common Use Cases

### 1. Documentation

```bash
# Generate diagram for README
python agent.py --save-graph ./docs/architecture/workflow

# Add to README.md
echo "![Workflow](./docs/architecture/workflow.png)" >> README.md
```

### 2. Debugging

```python
# In development/debugging script
from my_agent import create_graph, save_graph_visualization

graph = create_graph().compile()
save_graph_visualization(graph, f"./debug/workflow_v{version}")
```

### 3. Multiple Workflow Variants

```bash
# Save different configurations
python agent.py --mode hitl --save-graph ./graphs/hitl_mode
python agent.py --mode auto --save-graph ./graphs/auto_mode
```

### 4. CI/CD Pipeline

```yaml
# .github/workflows/docs.yml
- name: Generate workflow diagrams
  run: |
    python agent.py --save-graph ./docs/diagrams/workflow
    git add docs/diagrams/
    git commit -m "Update workflow diagrams"
```

---

## Troubleshooting

### Issue: "No module named 'pygraphviz'"

**Solution:** Install graphviz backend (optional, Mermaid works without it)
```bash
# macOS
brew install graphviz
pip install pygraphviz

# Ubuntu/Debian
sudo apt-get install graphviz graphviz-dev
pip install pygraphviz
```

### Issue: PNG is corrupted or empty

**Solution:** Check graph compilation
```python
# Verify graph compiles without errors
try:
    graph = builder.compile()
    print("Graph compiled successfully")
except Exception as e:
    print(f"Compilation error: {e}")
```

### Issue: Mermaid syntax errors when viewing

**Solution:** LangGraph generates valid Mermaid v10+ syntax. Update your viewer:
- Mermaid Live: Use https://mermaid.live (latest version)
- VS Code: Update Mermaid extension
- GitHub: Should work automatically

### Issue: Graph is too large/complex to read

**Solution:** Use `xray=False` for simplified view
```python
# Simplified view
png_data = graph.get_graph(xray=False).draw_mermaid_png()

# Or filter by node type in Mermaid editor after export
```

---

## Best Practices

1. **Version control:** Save graphs with version/date in filename
   ```bash
   python agent.py --save-graph ./docs/workflow_v2.1_20251224
   ```

2. **Organized storage:** Keep graphs in dedicated directory
   ```
   project/
   ├── docs/
   │   ├── diagrams/
   │   │   ├── workflow_main.png
   │   │   ├── workflow_main.mermaid
   │   │   └── workflow_debug.png
   ```

3. **Update regularly:** Regenerate after significant workflow changes

4. **Include in documentation:** Link or embed in README, wiki, or docs

5. **Compare versions:** Use git to track diagram changes over time

---

## Quick Reference

| Task | Command |
|------|---------|
| Save graph (CLI) | `python agent.py --save-graph name` |
| Save graph (Python) | `save_graph_visualization(graph, "name")` |
| Display in notebook | `display(Image(graph.get_graph(xray=True).draw_mermaid_png()))` |
| View PNG | `open workflow.png` |
| View Mermaid online | https://mermaid.live |
| Simplified graph | `graph.get_graph(xray=False)` |
| Detailed graph | `graph.get_graph(xray=True)` |

---

## Example: Complete Workflow

```python
#!/usr/bin/env python
"""Example: Agent with graph saving capability."""

from langgraph.graph import StateGraph, START, END
from typing import TypedDict
import argparse


class State(TypedDict):
    messages: list


def my_node(state: State):
    return {"messages": state["messages"] + ["processed"]}


def save_graph_visualization(graph, output_path: str = "workflow"):
    import os
    output_dir = os.path.dirname(output_path) or "."
    os.makedirs(output_dir, exist_ok=True)

    png_data = graph.get_graph(xray=True).draw_mermaid_png()
    with open(f"{output_path}.png", "wb") as f:
        f.write(png_data)

    mermaid_code = graph.get_graph(xray=True).draw_mermaid()
    with open(f"{output_path}.mermaid", "w") as f:
        f.write(mermaid_code)

    print(f"Saved: {output_path}.png and {output_path}.mermaid")


def create_workflow():
    builder = StateGraph(State)
    builder.add_node("process", my_node)
    builder.add_edge(START, "process")
    builder.add_edge("process", END)
    return builder


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--save-graph", type=str, help="Save graph to path")
    args = parser.parse_args()

    builder = create_workflow()
    graph = builder.compile()

    if args.save_graph:
        save_graph_visualization(graph, args.save_graph)
    else:
        # Run workflow
        result = graph.invoke({"messages": []})
        print(result)
```

**Usage:**
```bash
# Save graph
python example.py --save-graph ./my_workflow

# Run workflow
python example.py
```

---

## Related Resources

- [LangGraph Documentation](https://langchain-ai.github.io/langgraph/)
- [Mermaid Documentation](https://mermaid.js.org/)
- [Mermaid Live Editor](https://mermaid.live/)

---

**Note:** This guide is applicable to any LangGraph workflow, not just the HITL implementation.
