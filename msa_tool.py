import os
from itertools import combinations
from collections import Counter


class SequenceAlignmentTool:
    def __init__(self):
        self.protein_chars = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")
        self.nucleotide_chars = set("ATGCUN")

    # ---------------- Sequence Utilities ---------------- #
    def clean_sequence(self, seq):
        return seq.strip().replace(" ", "").replace("\n", "").upper()

    def validate_sequence(self, seq):
        seq = self.clean_sequence(seq)
        if not seq:
            raise ValueError("Sequence is empty.")
        filtered = seq.replace("*", "").replace("-", "")
        if not filtered.isalpha():
            raise ValueError(
                f"Sequence contains invalid characters: {seq[:30]}...")
        return seq.replace("-", "")

    def detect_sequence_type(self, seq):
        seq = self.clean_sequence(seq).replace("-", "").replace("*", "")
        chars = set(seq)
        if chars.issubset(self.nucleotide_chars):
            return "Nucleotide"
        return "Protein"

    def parse_fasta_text(self, text):
        sequences = []
        current_name = None
        current_seq = []

        for line in text.strip().splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    sequences.append((current_name, "".join(current_seq)))
                current_name = line[1:].strip() or f"Seq{len(sequences)+1}"
                current_seq = []
            else:
                current_seq.append(line)

        if current_name is not None:
            sequences.append((current_name, "".join(current_seq)))

        if not sequences:
            raise ValueError("No FASTA sequences found.")

        return [(name, self.validate_sequence(seq)) for name, seq in sequences]

    def load_fasta_file(self, filepath):
        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        with open(filepath, "r", encoding="utf-8") as f:
            return self.parse_fasta_text(f.read())

    # ---------------- Scoring ---------------- #
    def score_pair(self, a, b, match=2, mismatch=-1):
        if a == "-" or b == "-":
            return 0
        return match if a == b else mismatch

    # ---------------- Pairwise: Needleman-Wunsch (Global) with AFFINE gap ---------------- #
    def needleman_wunsch(self, seq1, seq2, match=2, mismatch=-1,
                         gap_open=-10, gap_extend=-0.5):
        """
        Global alignment using Needleman-Wunsch with AFFINE gap penalty.
        gap_open   : penalty to OPEN a new gap (paid once per gap run)
        gap_extend : penalty per residue EXTENDED in a gap (paid every position)

        Total gap cost for a gap of length k = gap_open + gap_extend * k
        This matches the EMBOSS Needle convention.
        """
        n, m = len(seq1), len(seq2)
        NEG_INF = float("-inf")

        # Three matrices: M (match/mismatch), X (gap in seq2), Y (gap in seq1)
        M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
        X = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
        Y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

        M[0][0] = 0.0
        for i in range(1, n + 1):
            X[i][0] = gap_open + gap_extend * i
        for j in range(1, m + 1):
            Y[0][j] = gap_open + gap_extend * j

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                s = match if seq1[i-1] == seq2[j-1] else mismatch

                best_prev = max(
                    M[i-1][j-1] if M[i-1][j-1] != NEG_INF else NEG_INF,
                    X[i-1][j-1] if X[i-1][j-1] != NEG_INF else NEG_INF,
                    Y[i-1][j-1] if Y[i-1][j-1] != NEG_INF else NEG_INF,
                )
                M[i][j] = (best_prev + s) if best_prev != NEG_INF else NEG_INF

                # Gap in seq2 (extend down in seq1)
                X[i][j] = max(
                    M[i-1][j] + gap_open +
                    gap_extend if M[i-1][j] != NEG_INF else NEG_INF,
                    X[i-1][j] + gap_extend if X[i-1][j] != NEG_INF else NEG_INF,
                    Y[i-1][j] + gap_open +
                    gap_extend if Y[i-1][j] != NEG_INF else NEG_INF,
                )

                # Gap in seq1 (extend right in seq2)
                Y[i][j] = max(
                    M[i][j-1] + gap_open +
                    gap_extend if M[i][j-1] != NEG_INF else NEG_INF,
                    X[i][j-1] + gap_open +
                    gap_extend if X[i][j-1] != NEG_INF else NEG_INF,
                    Y[i][j-1] + gap_extend if Y[i][j-1] != NEG_INF else NEG_INF,
                )

        # Final score: best of three matrices at (n, m)
        final_score = max(
            M[n][m] if M[n][m] != NEG_INF else NEG_INF,
            X[n][m] if X[n][m] != NEG_INF else NEG_INF,
            Y[n][m] if Y[n][m] != NEG_INF else NEG_INF,
        )

        # Traceback
        aligned1, aligned2 = [], []
        i, j = n, m

        # Determine starting matrix
        scores_at_end = {
            "M": M[n][m], "X": X[n][m], "Y": Y[n][m]
        }
        state = max(scores_at_end, key=scores_at_end.get)

        while i > 0 or j > 0:
            if state == "M":
                aligned1.append(seq1[i-1])
                aligned2.append(seq2[j-1])
                s = match if seq1[i-1] == seq2[j-1] else mismatch
                prev = M[i][j] - s
                if i > 0 and j > 0:
                    if abs((M[i-1][j-1] if M[i-1][j-1] != NEG_INF else NEG_INF) - prev) < 1e-9:
                        state = "M"
                    elif abs((X[i-1][j-1] if X[i-1][j-1] != NEG_INF else NEG_INF) - prev) < 1e-9:
                        state = "X"
                    else:
                        state = "Y"
                i -= 1
                j -= 1
            elif state == "X":
                aligned1.append(seq1[i-1])
                aligned2.append("-")
                if i > 0:
                    ext = X[i][j] - gap_extend
                    opn = X[i][j] - gap_open - gap_extend
                    if M[i-1][j] != NEG_INF and abs(M[i-1][j] - opn) < 1e-9:
                        state = "M"
                    elif X[i-1][j] != NEG_INF and abs(X[i-1][j] - ext) < 1e-9:
                        state = "X"
                    else:
                        state = "Y"
                i -= 1
            else:  # Y
                aligned1.append("-")
                aligned2.append(seq2[j-1])
                if j > 0:
                    ext = Y[i][j] - gap_extend
                    opn = Y[i][j] - gap_open - gap_extend
                    if M[i][j-1] != NEG_INF and abs(M[i][j-1] - opn) < 1e-9:
                        state = "M"
                    elif X[i][j-1] != NEG_INF and abs(X[i][j-1] - opn) < 1e-9:
                        state = "X"
                    else:
                        state = "Y"
                j -= 1

        a1 = "".join(reversed(aligned1))
        a2 = "".join(reversed(aligned2))
        return a1, a2, round(final_score, 1)

    # ---------------- Pairwise: Smith-Waterman (Local) with AFFINE gap ----------- #
    def smith_waterman(self, seq1, seq2, match=2, mismatch=-1,
                       gap_open=-10, gap_extend=-0.5):
        """
        Local alignment using Smith-Waterman with AFFINE gap penalty.
        Matches EMBOSS Water convention.
        """
        n, m = len(seq1), len(seq2)
        NEG_INF = float("-inf")

        M = [[0.0] * (m + 1) for _ in range(n + 1)]
        X = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
        Y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

        max_score = 0.0
        max_pos = (0, 0)

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                s = match if seq1[i-1] == seq2[j-1] else mismatch

                m_prev = max(M[i-1][j-1], X[i-1][j-1], Y[i-1][j-1])
                M[i][j] = max(0.0, m_prev + s)

                X[i][j] = max(
                    M[i-1][j] + gap_open + gap_extend,
                    X[i-1][j] + gap_extend,
                    Y[i-1][j] + gap_open + gap_extend,
                    0.0,
                )

                Y[i][j] = max(
                    M[i][j-1] + gap_open + gap_extend,
                    X[i][j-1] + gap_open + gap_extend,
                    Y[i][j-1] + gap_extend,
                    0.0,
                )

                best = max(M[i][j], X[i][j], Y[i][j])
                if best > max_score:
                    max_score = best
                    max_pos = (i, j)

        # Traceback
        aligned1, aligned2 = [], []
        i, j = max_pos

        scores_here = {"M": M[i][j], "X": X[i][j], "Y": Y[i][j]}
        state = max(scores_here, key=scores_here.get)

        while i > 0 and j > 0:
            cur = max(M[i][j], X[i][j], Y[i][j])
            if cur <= 0:
                break

            if state == "M":
                aligned1.append(seq1[i-1])
                aligned2.append(seq2[j-1])
                s = match if seq1[i-1] == seq2[j-1] else mismatch
                prev = M[i][j] - s
                p_m = M[i-1][j-1]
                p_x = X[i-1][j-1]
                p_y = Y[i-1][j-1]
                if abs(p_m - prev) < 1e-9:
                    state = "M"
                elif abs(p_x - prev) < 1e-9:
                    state = "X"
                else:
                    state = "Y"
                i -= 1
                j -= 1
            elif state == "X":
                aligned1.append(seq1[i-1])
                aligned2.append("-")
                ext = X[i][j] - gap_extend
                if abs(M[i-1][j] - ext + gap_open) < 1e-9:
                    state = "M"
                elif abs(X[i-1][j] - ext) < 1e-9:
                    state = "X"
                else:
                    state = "Y"
                i -= 1
            else:
                aligned1.append("-")
                aligned2.append(seq2[j-1])
                ext = Y[i][j] - gap_extend
                if abs(M[i][j-1] - ext + gap_open) < 1e-9:
                    state = "M"
                elif abs(X[i][j-1] - ext + gap_open) < 1e-9:
                    state = "X"
                else:
                    state = "Y"
                j -= 1

        a1 = "".join(reversed(aligned1))
        a2 = "".join(reversed(aligned2))
        return a1, a2, round(max_score, 1)

    # ---------------- Alignment Statistics ---------------- #
    def alignment_symbols(self, aln1, aln2):
        symbols = []
        matches = mismatches = gaps = 0
        for a, b in zip(aln1, aln2):
            if a == "-" or b == "-":
                symbols.append(" ")
                gaps += 1
            elif a == b:
                symbols.append("|")
                matches += 1
            else:
                symbols.append(".")
                mismatches += 1
        return "".join(symbols), matches, mismatches, gaps

    def compute_stats(self, aln1, aln2):
        """
        Stats matching EMBOSS convention:
          Identity  = matches / alignment_length * 100
          Similarity = (matches + similar) / alignment_length * 100
          Gaps       = gap_columns / alignment_length * 100
        """
        symbols, matches, mismatches, gaps = self.alignment_symbols(aln1, aln2)
        length = len(aln1)
        identity = (matches / length * 100) if length > 0 else 0.0
        similarity = ((matches + mismatches) / length *
                      100) if length > 0 else 0.0
        gap_pct = (gaps / length * 100) if length > 0 else 0.0
        return {
            "symbols": symbols,
            "matches": matches,
            "mismatches": mismatches,
            "gaps": gaps,
            "identity": identity,
            "similarity": similarity,
            "gap_percent": gap_pct,
            "length": length,
        }

    # ---------------- Multiple Sequence Alignment (Progressive) ---------------- #
    def pad_sequences(self, aligned_sequences):
        max_len = max(len(seq) for _, seq in aligned_sequences)
        return [(name, seq + "-" * (max_len - len(seq))) for name, seq in aligned_sequences]

    def _merge_aligned_with_profile(self, profile_seqs, new_name, new_seq):
        """
        Align new_seq against the consensus of profile_seqs, then
        project gaps back into all profile sequences.
        """
        max_len = max(len(s) for _, s in profile_seqs)
        consensus = []
        for col in range(max_len):
            col_chars = [s[col]
                         for _, s in profile_seqs if col < len(s) and s[col] != "-"]
            if col_chars:
                consensus.append(Counter(col_chars).most_common(1)[0][0])
            else:
                consensus.append("-")
        consensus_str = "".join(c for c in consensus if c != "-")

        aln_consensus, aln_new, _ = self.needleman_wunsch(
            consensus_str, new_seq)

        new_profile = list(profile_seqs)
        cons_idx = 0
        insert_positions = []

        for i, c in enumerate(aln_consensus):
            if c == "-":
                insert_positions.append(cons_idx)
            else:
                cons_idx += 1

        for pos in sorted(set(insert_positions)):
            new_profile = [(name, seq[:pos] + "-" + seq[pos:])
                           for name, seq in new_profile]

        new_profile.append((new_name, aln_new))
        new_profile = self.pad_sequences(new_profile)
        return new_profile

    def progressive_msa(self, sequences):
        """
        Progressive MSA using guide-tree ordering by pairwise similarity.
        Returns (aligned_sequences, pairwise_scores_dict).
        """
        if len(sequences) < 2:
            raise ValueError("At least two sequences are required for MSA.")

        validated = [(name, self.validate_sequence(seq))
                     for name, seq in sequences]

        n = len(validated)
        score_matrix = {}
        for i, j in combinations(range(n), 2):
            _, _, score = self.needleman_wunsch(
                validated[i][1], validated[j][1])
            score_matrix[(i, j)] = score

        avg_scores = []
        for i in range(n):
            pairs = [score_matrix.get((min(i, j), max(i, j)), 0)
                     for j in range(n) if j != i]
            avg_scores.append((sum(pairs) / len(pairs) if pairs else 0, i))
        order = [i for _, i in sorted(avg_scores, reverse=True)]

        profile = [validated[order[0]]]
        for idx in order[1:]:
            name, seq = validated[idx]
            profile = self._merge_aligned_with_profile(profile, name, seq)

        name_to_aligned = {name: seq for name, seq in profile}
        result = [(name, name_to_aligned[name]) for name, _ in validated]

        return result, score_matrix

    def msa_conservation(self, aligned_seqs):
        """Return per-column conservation symbol (ClustalW style)."""
        if not aligned_seqs:
            return ""
        length = max(len(s) for _, s in aligned_seqs)
        conservation = []
        for col in range(length):
            chars = set()
            has_gap = False
            for _, seq in aligned_seqs:
                if col < len(seq):
                    c = seq[col]
                    if c == "-":
                        has_gap = True
                    else:
                        chars.add(c)
            if has_gap or len(chars) == 0:
                conservation.append(" ")
            elif len(chars) == 1:
                conservation.append("*")
            else:
                conservation.append(".")
        return "".join(conservation)
