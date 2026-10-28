"""Regenerate local vector assets for the COSI README (Python standard library)."""
from pathlib import Path
from html import escape

out = Path(__file__).resolve().parent / "assets"
out.mkdir(parents=True, exist_ok=True)

banner = '''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="190" viewBox="0 0 1200 190" role="img" aria-labelledby="title desc">
<title id="title">Split data, estimate, combine</title>
<desc id="desc">A blue and teal header with a schematic of several interval estimates. Decorative intervals have no numerical scale.</desc>
<rect x=".5" y=".5" width="1199" height="189" rx="18" fill="#EDF6FA" stroke="#D6E6ED"/>
<rect x="820" y="20" width="350" height="150" rx="14" fill="#E2F3EC"/>
<path d="M18 2H450" stroke="#0072B2" stroke-width="4"/><path d="M450 2H770" stroke="#009E73" stroke-width="4"/><path d="M770 2H1020" stroke="#7B2CBF" stroke-width="4"/><path d="M1020 2H1182" stroke="#D55E00" stroke-width="4"/>
<g font-family="sans-serif">
 <text x="40" y="53" font-size="13" letter-spacing="2.5" font-weight="700" fill="#006F68">CONFIDENCE THROUGH STRUCTURE</text>
 <text x="40" y="99" font-size="32" font-weight="700" fill="#163A54">Split data. Estimate. Combine.</text>
 <text x="40" y="140" font-size="19" fill="#486578">Scalar intervals, vector sets, and monotone bands</text>
</g>
<g transform="translate(850 32)" fill="none" stroke-linecap="round">
 <path d="M124 0V124" stroke="#A8B9CA" stroke-dasharray="3 6" stroke-width="1.6"/>
 <g stroke="#0072B2" stroke-width="4"><path d="M35 14H210"/><path d="M35 7V21M210 7V21"/></g><circle cx="116" cy="14" r="6" fill="#0072B2"/>
 <g stroke="#009E73" stroke-width="4"><path d="M67 45H240"/><path d="M67 38V52M240 38V52"/></g><circle cx="145" cy="45" r="6" fill="#009E73"/>
 <g stroke="#7B2CBF" stroke-width="4"><path d="M48 76H188"/><path d="M48 69V83M188 69V83"/></g><circle cx="119" cy="76" r="6" fill="#7B2CBF"/>
 <g stroke="#D55E00" stroke-width="4"><path d="M80 107H222"/><path d="M80 100V114M222 100V114"/></g><circle cx="150" cy="107" r="6" fill="#D55E00"/>
</g>
</svg>'''
(out / "cosi_banner.svg").write_text(banner)

badges = [
    ("means", "Means", 94, "#005B8F", "#E8F3FC", "#B9DAEF"),
    ("regression", "Regression", 132, "#006F55", "#E8F7F0", "#B4DDCE"),
    ("dependence", "Dependence", 141, "#985100", "#FFF3E5", "#EBD0AF"),
    ("highdim", "High dimension", 161, "#6940A5", "#F2EBFC", "#D5C3ED"),
    ("monotone", "Monotone bands", 169, "#006C76", "#E6F5F7", "#B8DBDF"),
]
for name, label, width, ink, background, border in badges:
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="28" viewBox="0 0 {width} 28" role="img" aria-label="{escape(label)}">
<rect x=".5" y=".5" width="{width-1}" height="27" rx="7" fill="{background}" stroke="{border}"/>
<circle cx="14" cy="14" r="3.8" fill="{ink}"/>
<text x="25" y="18.5" fill="{ink}" font-family="DejaVu Sans, Arial, sans-serif" font-size="12" font-weight="600">{escape(label)}</text>
</svg>'''
    (out / f"badge_{name}.svg").write_text(svg)
print("Created banner and five linked-feature badge assets.")
