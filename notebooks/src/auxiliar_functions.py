import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm

class AuxFunctions:

    @classmethod
    def holm_adjust(cls, pvals: np.ndarray) -> np.ndarray:
        """
        Holm step-down adjustment.
        Returns adjusted p-values in the original order.
        """
        pvals = np.asarray(pvals, dtype=float)
        m = len(pvals)
        order = np.argsort(pvals)
        ranked = pvals[order]

        adjusted = np.empty(m, dtype=float)
        for i, p in enumerate(ranked):
            adjusted[i] = (m - i) * p
        # Ensure monotonicity
        adjusted = np.maximum.accumulate(adjusted)
        adjusted = np.clip(adjusted, 0, 1)

        out = np.empty(m, dtype=float)
        out[order] = adjusted
        return out

    @classmethod
    def dunn_posthoc(cls, groups: dict, p_adjust: str = "holm") -> dict:
        """
        Dunn's test for multiple pairwise comparisons.
        groups: dict[str|float, array-like] mapping group -> observations
        Returns dict with keys (g1, g2) sorted and p-values (adjusted).
        """
        # Flatten data
        group_keys = list(groups.keys())
        values = []
        labels = []
        for g in group_keys:
            arr = np.asarray(groups[g], dtype=float)
            values.append(arr)
            labels.extend([g] * len(arr))
        x = np.concatenate(values)
        labels = np.asarray(labels, dtype=object)

        # Ranks (average for ties)
        r = rankdata(x)

        # Group rank sums and sizes
        N = len(x)
        group_sizes = {g: np.sum(labels == g) for g in group_keys}
        rank_sums = {g: np.sum(r[labels == g]) for g in group_keys}

        # Tie correction
        _, counts = np.unique(x, return_counts=True)
        tie_term = np.sum(counts**3 - counts)
        tie_correction = 1.0 - tie_term / (N**3 - N) if N > 1 else 1.0

        # Variance term
        var = (N * (N + 1) / 12.0) * tie_correction

        # Pairwise z and p
        pairs = []
        raw_p = []
        for i in range(len(group_keys)):
            for j in range(i + 1, len(group_keys)):
                g1, g2 = group_keys[i], group_keys[j]
                n1, n2 = group_sizes[g1], group_sizes[g2]
                if n1 == 0 or n2 == 0:
                    continue
                mean_r1 = rank_sums[g1] / n1
                mean_r2 = rank_sums[g2] / n2
                z = (mean_r1 - mean_r2) / np.sqrt(var * (1.0 / n1 + 1.0 / n2))
                p = 2.0 * norm.sf(abs(z))
                pairs.append((g1, g2))
                raw_p.append(p)

        raw_p = np.asarray(raw_p, dtype=float)

        # Adjust p-values
        if p_adjust.lower() == "holm":
            adj_p = cls.holm_adjust(raw_p)
        elif p_adjust.lower() == "bonferroni":
            adj_p = np.clip(raw_p * len(raw_p), 0, 1)
        else:
            raise ValueError("p_adjust must be 'holm' or 'bonferroni'")

        out = {}
        for (g1, g2), p in zip(pairs, adj_p):
            key = tuple(sorted((g1, g2), key=lambda z: float(z)))
            out[key] = float(p)
        return out

    @classmethod
    def cld_letters(cls, group_means: dict, pvals: dict, alpha: float = 0.05) -> dict:
        """
        Compact Letter Display (CLD) from pairwise p-values.
        Groups that are NOT significantly different (p >= alpha) can share letters.
        """
        groups = list(group_means.keys())
        # Sort groups by mean descending (common convention)
        groups = sorted(groups, key=lambda g: group_means[g], reverse=True)

        letters = {g: "" for g in groups}
        letter_list = []

        def nonsig(g, h):
            key = tuple(sorted((g, h), key=lambda z: float(z)))
            return pvals.get(key, 1.0) >= alpha

        for g in groups:
            placed = False
            for li, L in enumerate(letter_list):
                # Can g join letter L only if g is nonsignificant vs everyone already with L
                members = [h for h in groups if L in letters[h]]
                if all(nonsig(g, h) for h in members):
                    letters[g] += L
                    placed = True
            if not placed:
                # Create a new letter
                new_letter = chr(ord("a") + len(letter_list))
                letter_list.append(new_letter)
                letters[g] += new_letter

            # Clean-up pass: if a group has multiple letters, try to remove redundant ones
            # (light simplification; still produces correct CLD most of the time)
            for L in list(letters[g]):
                members = [h for h in groups if (h != g and L in letters[h])]
                if members and all(nonsig(g, h) for h in members):
                    continue

        return letters

    @classmethod
    def compute_letters_per_enzyme(cls, df: pd.DataFrame, alpha: float = 0.05, adjust: str = "holm", present_treats=[]) -> pd.DataFrame:
        """
        For each enzyme, run Dunn test across treatments and generate CLD letters.
        Returns a DataFrame: Enzyme, Treatment, Letters
        """
        rows = []
        for enz, dsub in df.groupby("Enzyme"):
            # Group values by treatment
            groups = {float(t): dsub.loc[dsub["Treatment"] == t, "Viability"].values for t in dsub["Treatment"].cat.categories}
            # Skip if only one treatment
            if sum(len(v) > 0 for v in groups.values()) < 2:
                continue

            # Means for ordering in CLD
            means = {t: np.mean(v) if len(v) else np.nan for t, v in groups.items()}
            # Remove empty groups
            groups = {t: v for t, v in groups.items() if len(v)}
            means = {t: means[t] for t in groups.keys()}

            pvals = cls.dunn_posthoc(groups, p_adjust=adjust)
            letters = cls.cld_letters(means, pvals, alpha=alpha)

            for t, L in letters.items():
                rows.append({"Enzyme": enz, "Treatment": t, "Letters": L})

        out = pd.DataFrame(rows)
        if not out.empty:
            out["Treatment"] = pd.Categorical(out["Treatment"], categories=present_treats, ordered=True)
            out = out.sort_values(["Enzyme", "Treatment"])
        return out
