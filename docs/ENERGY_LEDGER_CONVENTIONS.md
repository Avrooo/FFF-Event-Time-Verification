# Energy-ledger conventions and reference states

The archival `N0/N0_results/A0_md0p1_audit.json` reports:

- `heat_in_J` = 1463.0795379612719 J (newly deposited thermal excess above ambient)
- `stored_change_J` = -1247.3802958205176 J (change from newly deposited states to the final states)
- `natural_conv_J` = 393.26075375749423 J (outward)
- `radiation_J` = 987.8732777122755 J (outward)
- `bed_net_J` = -120.57432425870198 J (**outward positive** in the archive, i.e., negative means bed heats the part)
- `balance_error_J` = +13.179411390550058 J (archive convention)
- absolute discrepancy = 0.9007993788850952%.

The manuscript uses **bed exchange positive INTO the part**, thus

`Qbed_manuscript = -bed_net_J = +120.57432425870198 J`.

The manuscript's stored-energy variable is final thermal excess above ambient, not the archived change-from-deposition value:

`DeltaU_manuscript = heat_in_J + stored_change_J = 215.6992421407543 J`.

With the manuscript definition

`R_manuscript = heat_in_J + Qbed_manuscript - natural_conv_J - radiation_J - DeltaU_manuscript`

we get `R_manuscript = -13.179411390550058 J`, the negative of the archived `balance_error_J` by definition. Hence the same absolute discrepancy `100*abs(R_manuscript)/heat_in_J = 0.9007993788850952%`.

**Do not flip signs in the raw JSON/CSV files**, which retain their original conventions. When presenting values or regenerating the paper, apply this documented conversion consistently. The archived 0.9008% refers to reconstructed physical-clock accounting, not a comparison with printer measurements.
