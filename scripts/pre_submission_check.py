#!/usr/bin/env python3
"""
Pre-submission validation script for Nature Genetics manuscript.

Checks:
1. Citation key consistency between main.tex and references.bib
2. Figure references and includes
3. Required sections (Data availability, Code availability)
4. No placeholder text remaining
5. Compiles without errors (optional)

Usage:
    python scripts/pre_submission_check.py
    python scripts/pre_submission_check.py --compile  # Also run latexmk
"""

import re
import subprocess
import sys
from pathlib import Path

# Paths
MANUSCRIPT_DIR = Path(__file__).parent.parent / "manuscript"
MAIN_TEX = MANUSCRIPT_DIR / "main.tex"
REFERENCES_BIB = MANUSCRIPT_DIR / "references.bib"


def extract_citations_from_tex(tex_path: Path) -> set[str]:
    """Extract all citation keys used in the tex file."""
    content = tex_path.read_text(encoding="utf-8")
    # Match \cite{key1,key2,...} patterns
    pattern = r"\\cite\{([^}]+)\}"
    citations = set()
    for match in re.finditer(pattern, content):
        keys = match.group(1).split(",")
        citations.update(k.strip() for k in keys)
    return citations


def extract_keys_from_bib(bib_path: Path) -> set[str]:
    """Extract all entry keys defined in the bib file."""
    content = bib_path.read_text(encoding="utf-8")
    # Match @type{key, patterns
    pattern = r"@\w+\{(\w+),"
    return set(re.findall(pattern, content))


def check_citation_consistency():
    """Check that all cited keys exist in the bibliography."""
    print("\n" + "=" * 60)
    print("CHECKING CITATION CONSISTENCY")
    print("=" * 60)
    
    tex_citations = extract_citations_from_tex(MAIN_TEX)
    bib_keys = extract_keys_from_bib(REFERENCES_BIB)
    
    missing = tex_citations - bib_keys
    unused = bib_keys - tex_citations
    
    if missing:
        print(f"\n❌ MISSING in references.bib (cited but not defined):")
        for key in sorted(missing):
            print(f"   - {key}")
        return False
    else:
        print(f"\n✓ All {len(tex_citations)} cited keys found in references.bib")
    
    if unused:
        print(f"\n⚠ UNUSED in references.bib (defined but not cited):")
        for key in sorted(unused):
            print(f"   - {key}")
    
    return len(missing) == 0


def check_figure_references():
    """Check figure labels and references are consistent."""
    print("\n" + "=" * 60)
    print("CHECKING FIGURE REFERENCES")
    print("=" * 60)
    
    content = MAIN_TEX.read_text(encoding="utf-8")
    
    # Extract figure labels
    labels = set(re.findall(r"\\label\{(fig:[^}]+)\}", content))
    # Extract figure references
    refs = set(re.findall(r"\\ref\{(fig:[^}]+)\}", content))
    
    missing_labels = refs - labels
    unreferenced = labels - refs
    
    if missing_labels:
        print(f"\n❌ UNDEFINED figure labels (referenced but not defined):")
        for label in sorted(missing_labels):
            print(f"   - {label}")
        return False
    else:
        print(f"\n✓ All {len(refs)} figure references have corresponding labels")
    
    if unreferenced:
        print(f"\n⚠ UNREFERENCED figures (defined but not referenced):")
        for label in sorted(unreferenced):
            print(f"   - {label}")
    
    # Check for includegraphics (even if commented)
    includes = re.findall(r"%?\s*\\includegraphics.*?\{([^}]+)\}", content)
    print(f"\n✓ Found {len(includes)} figure include statements")
    
    return len(missing_labels) == 0


def check_required_sections():
    """Check that required Nature sections are present."""
    print("\n" + "=" * 60)
    print("CHECKING REQUIRED SECTIONS")
    print("=" * 60)
    
    content = MAIN_TEX.read_text(encoding="utf-8")
    
    required = [
        (r"\\section\*?\{Data availability\}", "Data availability"),
        (r"\\section\*?\{Code availability\}", "Code availability"),
        (r"\\section\{Methods\}", "Methods"),
        (r"\\section\{Results\}", "Results"),
        (r"\\section\{Discussion\}", "Discussion"),
        (r"\\begin\{abstract\}", "Abstract"),
    ]
    
    all_present = True
    for pattern, name in required:
        if re.search(pattern, content, re.IGNORECASE):
            print(f"✓ {name} section found")
        else:
            print(f"❌ {name} section MISSING")
            all_present = False
    
    return all_present


def check_placeholders():
    """Check for remaining placeholder text."""
    print("\n" + "=" * 60)
    print("CHECKING FOR PLACEHOLDERS")
    print("=" * 60)
    
    content = MAIN_TEX.read_text(encoding="utf-8")
    
    # Common placeholder patterns
    placeholders = [
        (r"\[TODO[^\]]*\]", "TODO markers"),
        (r"\[INSERT[^\]]*\]", "INSERT markers"),
        (r"\[PLACEHOLDER[^\]]*\]", "PLACEHOLDER markers"),
        (r"XXX", "XXX markers"),
        (r"\?\?\?", "??? markers"),
        (r"\.{3,}", "... ellipsis (may be intentional)"),
    ]
    
    issues = []
    for pattern, name in placeholders:
        matches = re.findall(pattern, content, re.IGNORECASE)
        if matches:
            issues.append((name, len(matches)))
    
    if issues:
        print("\n⚠ Potential placeholders found:")
        for name, count in issues:
            print(f"   - {count}x {name}")
        return False
    else:
        print("\n✓ No obvious placeholder text found")
        return True


def check_compilation(compile_flag: bool):
    """Optionally check that the document compiles."""
    if not compile_flag:
        return True
    
    print("\n" + "=" * 60)
    print("CHECKING COMPILATION")
    print("=" * 60)
    
    try:
        result = subprocess.run(
            ["latexmk", "-pdf", "-interaction=nonstopmode", str(MAIN_TEX)],
            cwd=MANUSCRIPT_DIR,
            capture_output=True,
            text=True,
            timeout=120,
        )
        
        # Check for common errors
        log_file = MAIN_TEX.with_suffix(".log")
        if log_file.exists():
            log_content = log_file.read_text(encoding="utf-8", errors="ignore")
            
            # Check for undefined citations
            undefined_cites = re.findall(r"Citation `([^']+)' on page", log_content)
            if undefined_cites:
                print(f"\n❌ Undefined citations in compiled PDF:")
                for cite in set(undefined_cites):
                    print(f"   - {cite}")
                return False
            
            # Check for missing figures
            missing_figs = re.findall(r"File `([^']+)' not found", log_content)
            if missing_figs:
                print(f"\n⚠ Missing figure files (may be intentional):")
                for fig in set(missing_figs):
                    print(f"   - {fig}")
        
        if result.returncode == 0:
            print("\n✓ Document compiles successfully")
            return True
        else:
            print("\n❌ Compilation failed")
            print(result.stderr[:1000] if result.stderr else "No error output")
            return False
            
    except FileNotFoundError:
        print("\n⚠ latexmk not found - skipping compilation check")
        return True
    except subprocess.TimeoutExpired:
        print("\n⚠ Compilation timed out")
        return False


def main():
    """Run all pre-submission checks."""
    compile_flag = "--compile" in sys.argv
    
    print("\n" + "=" * 60)
    print("NATURE GENETICS PRE-SUBMISSION VALIDATION")
    print("=" * 60)
    print(f"Manuscript: {MAIN_TEX}")
    print(f"Bibliography: {REFERENCES_BIB}")
    
    results = [
        ("Citation consistency", check_citation_consistency()),
        ("Figure references", check_figure_references()),
        ("Required sections", check_required_sections()),
        ("Placeholder text", check_placeholders()),
    ]
    
    if compile_flag:
        results.append(("Compilation", check_compilation(True)))
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    all_passed = True
    for name, passed in results:
        status = "✓ PASS" if passed else "❌ FAIL"
        print(f"{status}: {name}")
        if not passed:
            all_passed = False
    
    if all_passed:
        print("\n🎉 All checks passed! Ready for submission.")
        return 0
    else:
        print("\n⚠ Some checks failed. Please review before submission.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
