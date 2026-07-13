# methods.md discipline

Maintain a running `methods.md` in the project root. For EVERY analytical step append:

```markdown
## <date> — <step name>
- **Method:** what was done (tool + version).
- **References:** PMID/DOI for the method.
- **Parameters & rationale:** key params and WHY these values.
- **Results:** interpretation of what was observed.
- **Troubleshooting:** for unexpected results, what was tried and the outcome.
```

This is the manuscript methods section accruing in real time. Never defer it.

## Self-documenting runs

Make runs reconstructable from disk, not just from memory:

- **Encode key params in output filenames** for swept steps, e.g.
  `com.hvg2000_PC10_res0.6.rds` or `model_lr1e-3_bs256_ep50.pt`. The filename then states the run.
- **Write a per-unit run/QC log table** (one row per sample / run / fold: inputs, key counts,
  thresholds applied, outputs) next to the results — turns a many-unit batch into an auditable
  table instead of buried console output.
