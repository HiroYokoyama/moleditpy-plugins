# Chat with Molecule (Gemini)

A plugin for MoleditPy that enables conversational interaction with the currently loaded molecule using Google's Gemini AI.

![Chat with Molecule 1](img/1.png)
![Chat with Molecule 2](img/2.png)
![Chat with Molecule 3](img/3.png)
![Chat with Molecule 4](img/4.png)

## Features

- **Context Awareness**: Automatically injects the SMILES string of the currently loaded molecule into the chat context.
- **Molecular Structure Recognition**: Instructions the AI to format molecular names as clickable links `[Name](smiles:...)` which automatically load the structure into the editor when clicked.
- **LaTeX Math Rendering**: Displays mathematical equations and chemical formulas using `matplotlib` (rendered as clear, scalable images).
- **Smart Link Handling**: Differentiates between:
    - **SMILES Links**: Load molecules directly in MoleditPy.
    - **Web Links**: Open in your default system browser (Chrome/Edge/etc.).
- **Model Selection**: Choose from available Gemini models (e.g., `gemini-1.5-flash`, `gemini-1.5-pro`).
- **History Export**: Export your chat session to a Markdown file.

## Installation

This plugin requires the `google-generativeai`, `markdown`, and `matplotlib` Python packages.

1.  Open your terminal or command prompt.
2.  Run the following command:

```bash
pip install google-generativeai markdown matplotlib
```

## Setup

1.  **Get an API Key**: Go to [Google AI Studio](https://aistudio.google.com/app/api-keys) and generate a free API key.
2.  **Enter Key**: Open the plugin in MoleditPy, paste your key into the "API Key" field, and click "Save & Reload".

## Usage

1.  Load a molecule in MoleditPy (or draw one).
2.  Open **Plugins > Chat with Molecule**.
3.  The "Context" label at the bottom will update to show the current SMILES string.
4.  Type your question (e.g., "What are the medicinal properties of this molecule?") and press Enter.
5.  **Equations**: Ask for formulas (e.g. "Show the reaction mechanism"). The plugin renders LaTeX safely (automatically simplifying complex tags for the viewer).
