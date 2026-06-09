# NA

------------------------------------------------------------------------

`ColFn()` in FarFirst:

Sometimes it’s simplest to compute self-distance; other times it’s
unnecessary.

Allow the user to decide. If the returned value is length N, assume that
the self distance has been reported, and ignore it. If the returned
value is length N - 1, assume the self distance has been omitted.

------------------------------------------------------------------------

In FarFirst(), when n \> N, run the algo anyway. The order of entries is
informative; and this allows a single search to be conducted and then
downsampled by selecting elements 1:.

Add a test such that FarFirst(n = N)\[1:5\] == FarFirst(n = 5).

------------------------------------------------------------------------
