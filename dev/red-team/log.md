# Red-team rotation log — MaxMin

Per-round notes for the `/red-team` skill. Newest entries appended at the bottom.
The `last_focus:` pointer at the very bottom tells the next invocation which area was last reviewed;
the next area is `(last_focus mod N) + 1` where N = number of areas in `focus-areas.md` (currently 7).

## Per-round entry format

```
### Round <n> — area #<N> <name>  (<YYYY-MM-DD>)
- tier: <sonnet|opus|fable>            # model the finder ran at
- finder yield: <count>                # confirmed non-trivial findings (after verification)
- verifier verdicts: <e.g. 3 candidates → 2 REAL, 1 REFUTED>
- trivial fixes: <count or list>
- seam: <"still yielding" | "ran dry">
- high-sev signal: <yes/no — detail if yes>
- next-visit decision: <"stay at <tier>" | "escalate to <tier>" | "escalate (high-sev signal)">
```

The skill reads the most recent entry for an area to set `last_tier:` and `yield:` for the escalation
decision. Escalation rules (SKILL.md §Model tiers / Normal run):
- yielded → re-visit SAME tier with a fresh agent (mine the seam)
- ran dry → escalate one tier (sonnet → opus → fable)
- high-sev signal → escalate immediately regardless

---

(no rounds yet)

---

last_focus: 0
