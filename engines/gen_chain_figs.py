from chain_diagrams import chain, txt

def kerr_chain():
    o = [txt(r"$\avg{\bd b}$", -3.7, 0), txt(r"$=$", -2.85, 0)]
    o.append(chain(["ld bd", "ld b"], [], ox=0.0, right=r"$r_{00}$"))
    o.append(txt(r"$+$", 1.4, 0))
    o.append(chain(["ld bd", "ld b"], [["ld bd", "ld b"]], ox=3.6,
                   right=r"$r_{00}$", seg_labels=[r"$O_{kl}$"]))
    o.append(txt(r"$+$", 8.8, 0))
    o.append(chain(["ld bd", "ld b"], [["ld bd", "ld b"], ["ld bd", "ld b"]], ox=11.0,
                   right=r"$r_{00}$", seg_labels=[r"$O_{kl}$", r"$O_{k'l'}$"]))
    o.append(txt(r"$+\;\cdots$", 18.4, 0))
    o.append(txt(r"$N=1$", 0.0, -1.3, "footnotesize"))
    o.append(txt(r"$N=2$", 4.75, -1.3, "footnotesize"))
    o.append(txt(r"$N=3$", 13.3, -1.3, "footnotesize"))
    o.append(txt(r"Fig.~\ref{fig:kerrdiagrams}(a)--(f)", 0.0, -1.85, "scriptsize"))
    o.append(txt(r"one factor per vertex, intermediate label summed", 9.0, -2.5, "scriptsize"))
    return "\n".join(o)

def ladder_chain():
    S = ["ld sp", "ld sm"]           # sigma^+ and sigma^- lines: the label Otilde_11
    C = ["ld sm"]                    # a single sigma^- line: the coherence Otilde_01
    o = [txt(r"$\avg{\sigp\sigm}$", -4.0, 0), txt(r"$=$", -3.0, 0)]
    o.append(chain(S, [C], ox=0.0, right=r"$\ketbra{g}{g}$", dx=2.1,
                   seg_labels=[r"$\tilde O_{01}$"]))
    o.append(txt(r"$+$", 3.9, 0))
    o.append(chain(S, [C, S, C], ox=6.0, right=r"$\ketbra{g}{g}$", dx=2.1,
                   seg_labels=[r"$\tilde O_{01}$", r"$\tilde O_{11}$", r"$\tilde O_{01}$"]))
    o.append(txt(r"$+\;\cdots\;+$", 14.4, 0))
    o.append(chain(S, [C, S, C, S, C], ox=17.2, right=r"$\ketbra{g}{g}$", dx=2.1,
                   dots_after=2,
                   seg_labels=[r"$\tilde O_{01}$", r"$\tilde O_{11}$", None,
                               r"$\tilde O_{11}$", r"$\tilde O_{01}$"]))
    o.append(txt(r"$+\;\cdots$", 30.4, 0))
    o.append(txt(r"$\dfrac{\Omega^2}{4D}$", 1.05, -1.45, "footnotesize"))
    o.append(txt(r"$\dfrac{\Omega^2}{4D}\Big(\!-\dfrac{\Omega^2}{2D}\Big)$", 9.2, -1.5, "footnotesize"))
    o.append(txt(r"$\dfrac{\Omega^2}{4D}\Big(\!-\dfrac{\Omega^2}{2D}\Big)^{\!r-1}$", 22.4, -1.5, "footnotesize"))
    o.append(txt(r"$N=2$", 1.05, -2.2, "scriptsize"))
    o.append(txt(r"$N=4$", 9.2, -2.2, "scriptsize"))
    o.append(txt(r"$N=2r$", 22.4, -2.2, "scriptsize"))
    return "\n".join(o)

open("tikz_kerrchain.tex","w").write(kerr_chain())
open("tikz_ladderchain.tex","w").write(ladder_chain())
print("ok")
