# Figure Reference Assets

These two PNG files are the high-resolution base panels used to reconstruct
Figure 1a and Figure 1b for manuscript figure revision 1.2.0-r5:

- `Figure1a_0806_base.png` (`2305 x 946` pixels)
- `Figure1b_0806_base.png` (`1266 x 791` pixels)

They preserve the accepted card artwork and layout. The public figure builder
replaces only versioned benchmark values before composing Figure 1. Keeping the
base panels here avoids a dependency on an unpublished manuscript DOCX and
makes `python paper/generate_manuscript_figures.py` runnable from a clean clone.
