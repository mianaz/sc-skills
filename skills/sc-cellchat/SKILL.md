---
name: sc-cellchat
description: Use when inferring ligand-receptor cell-cell communication from annotated scRNA-seq (CellChat).
---

# Single-Cell Cell–Cell Communication (CellChat v2)

## When to use
- User wants to know which cell types are "talking" to each other via secreted/membrane ligands and receptors
- User wants to compare communication networks across conditions (e.g., tumour vs normal)
- User wants pathway-level signalling summaries (e.g., MHC-II, VEGF, WNT signalling strength)
- User wants chord diagrams or bubble plots of interaction counts/weights

## When NOT to use
- Cell types are not yet annotated → run `sc-annotation` first; CellChat requires cell-type labels
- Fewer than ~200 cells per cell type → statistical power too low for reliable inference
- User wants physical proof of interaction → CellChat is statistical inference only; results are hypotheses

## Method decision tree
1. Is the Seurat/AnnData object annotated with cell-type labels? No → `sc-annotation` first.
2. Single condition → standard CellChat v2 workflow: createCellChat → filterCommunication → computeCommunicationProb → aggregateNet.
3. Two conditions to compare → run CellChat on each separately, then liftCellChat + mergeCellChat → compareInteractions.
4. Want pathway-level (not pairwise) summaries? → netAnalysis_computeCentrality + signalling pathway heatmap.

## Key caveats
- All inferences are based on co-expression of ligand in sender and receptor in receiver — statistical, not experimentally validated
- CellChat's built-in database (CellChatDB) is human/mouse only; verify species before running
- Treat chord diagram "interaction counts" as a starting point for hypothesis generation, not a conclusion
- High interaction count between two rare cell types may reflect noise from small cell numbers

## Routes to
- No direct scientific-agent-skill; CellChat v2 is R-only and self-contained
- For the upstream annotated object → `sc-annotation`
- For figure styling of chord/bubble plots → `scientific-plotting`
