import json, numpy as np
R = json.load(open('truncation_results.json')); t5 = json.load(open('table5_corrected.json'))
out = {}
out["kerr"] = {k: dict(N_star=v["N_star"], value=v["value"], delta=v["delta"], true=v["true_error"], N_true_opt=v["N_true_optimum"], true_opt=min(v["true_error_all_orders"]), exact=v["exact"]) for k, v in R["kerr_cavity"]["results"].items()}
rg = R["ring_K3"]; out["ring"] = {}
for key, r in rg["results"].items():
    e = rg["exact"][key]; r6 = rg["results_Nmax6"][key]
    out["ring"][key] = dict(g2=[r["g2"]["N_star"], r["g2"]["value"], r["g2"]["delta"], r["g2"]["true_error"]], g3=[r["g3"]["N_star"], r["g3"]["value"], r["g3"]["delta"], r["g3"]["true_error"]],
        g2_N6=[r6["g2"]["N_star"], r6["g2"]["delta"], r6["g2"]["true_error"]], g3_N6=[r6["g3"]["N_star"], r6["g3"]["delta"], r6["g3"]["true_error"]],
        g2_N2_err=r["g2"]["true_error_all_orders"][2], g2_N4_err=r["g2"]["true_error_all_orders"][4],
        exact=[e["g2"], e["g3"], e.get("g2_cutoff_uncertainty"), e.get("g3_cutoff_uncertainty"), e.get("source")], lab7=[e["lab_basis"]["g2"], e["lab_basis"]["g3"]], cached_g2=e.get("cached_ring_exactU", {}).get("g2"))
out["crosscheck"] = json.load(open('exact_pool_raw.json'))  # placeholder removed below
cc = [v for k, v, dt in out["crosscheck"] if k[0] == "series"]; out.pop("crosscheck")
k3 = t5["K3"]["0.05"]
out["K3"] = dict(exact=k3["exact"], exact_lab=k3["exact_lab_basis"], mf0=k3["mean_field_U0_displacement"], mfsc=k3["mean_field_selfconsistent"])
for k, v in k3.items():
    if k.startswith("diagrams"):
        out["K3"][k] = {q: dict(Nstar=v[q]["N_star"], value=v[q]["value"], delta=v[q]["delta"], true=v[q].get("true_error"), ps=[round(x, 6) for x in v[q]["partial_sums"]]) for q in ("c1", "fano", "c3c1")}
out["K1"] = {}
for U, d in t5["K1"].items():
    o = dict(exact=d["exact"], mf=d["mean_field"], gauss=d["gaussian"], gauss_contour=d["gaussian_contour"], disp=d.get("exact_displaced_frame"))
    for k, v in d.items():
        if k.startswith("diagrams_maxdeg14") and k.endswith("r0.1"):
            o[k] = {q: dict(Nstar=v[q]["N_star"], value=v[q]["value"], delta=v[q]["delta"], true=v[q]["true_error"], S3=v[q]["partial_sums"][3], S4=v[q]["partial_sums"][4]) for q in ("c1", "fano", "c3c1")}
    out["K1"][U] = o
fc = json.load(open('fig_fcs_errorbars.json'))
out["chain8"] = {k: dict(nb={q: [v["notebookRS"][q]["N_star"], v["notebookRS"][q]["value"], v["notebookRS"][q]["delta"]] for q in ("c1", "fano", "c3c1")},
                         rs={q: [v["fullRS"][q]["N_star"], v["fullRS"][q]["value"], v["fullRS"][q]["delta"], v["fullRS"][q]["partial_sums"][3]] for q in ("c1", "fano", "c3c1")},
                         gauss=v["gaussian"]["value"]) for k, v in fc["results"].items() if v.get("fullRS")}
tm = json.load(open('timing_results.json')); out["timing_fits"] = tm["fits"]
out["timing_pass3"] = {k: dict(K=v["K"], t=[round(x, 4) for x in v["t"]]) for k, v in tm["passes"]["pass3"].items()}
out["timing_pass1"] = {k: dict(K=v["K"], t=[round(x, 4) for x in v["t"]], load=[round(x) for x in v["load1"]]) for k, v in tm["passes"]["pass1_2"].items()}
print(json.dumps(out, indent=0, default=str))
