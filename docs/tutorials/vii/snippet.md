## Basis crushing
<div class="warning">
<b>WARNING:</b>
<br>
This feature is still under development and may not behave as expected.
<br>
</div>
<div>
<br>
</div>

You can "crush" the basis functions into a particular resolution so your analytical expression becomes shorter.
For example, on the $$12\times 12$$ grid for our example, we get 144 grid points, which results in many basis functions.
However, as we saw when analyzing the symmetry, many of those points are equivalent, so we should only fewer basis functions.
The `--crush` flag takes a tolerance that specifies how different the basis functions have to be to be considered an idependent function.
If basis functions are similar enough in value, their coefficients are combined, and they are reported as a single function.

