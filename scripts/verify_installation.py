#!/usr/bin/env python3
"""
Verify gbcms v2 installation and dependencies.

This script checks that all required components are properly installed
and can be imported.
"""

import importlib.util
import sys


def check_import(module_name: str, package_name: str | None = None) -> tuple[bool, str]:
    """
    Check if a module can be imported.

    Args:
        module_name: Name of the module to import
        package_name: Display name (defaults to module_name)

    Returns:
        Tuple of (success, message)
    """
    package_name = package_name or module_name
    try:
        if importlib.util.find_spec(module_name) is not None:
            return True, f"✅ {package_name}"
        else:
            return False, f"❌ {package_name}: Not found"
    except ImportError as e:
        return False, f"❌ {package_name}: {str(e)}"
    except Exception as e:
        return False, f"❌ {package_name}: Error {str(e)}"


def main():
    """Run installation verification."""
    print("=" * 60)
    print("gbcms v2 Installation Verification")
    print("=" * 60)
    print()

    # Core dependencies
    print("Core Dependencies:")
    print("-" * 60)

    core_deps = [
        ("pysam", "pysam"),
        ("typer", "typer"),
        ("rich", "rich"),
        ("pydantic", "pydantic"),
    ]

    core_results = [check_import(mod, pkg) for mod, pkg in core_deps]
    for _success, msg in core_results:
        print(msg)

    print()

    # gbcms modules
    print("gbcms Modules:")
    print("-" * 60)

    gb_modules = [
        ("gbcms", "gbcms package"),
        ("gbcms.cli", "CLI"),
        ("gbcms.models.core", "Core Models"),
        ("gbcms.io.input", "Input Readers"),
        ("gbcms.io.output", "Output Writers"),
        ("gbcms.pipeline", "Pipeline"),
        ("gbcms._rs", "Rust Extension (gbcms._rs)"),
    ]

    gb_results = [check_import(mod, pkg) for mod, pkg in gb_modules]
    for _success, msg in gb_results:
        print(msg)

    print()

    # Summary
    print("=" * 60)
    print("Summary:")
    print("-" * 60)

    all_results = core_results + gb_results
    total = len(all_results)
    passed = sum(1 for success, _ in all_results if success)
    failed = total - passed

    print(f"Total checks: {total}")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")

    if failed == 0:
        print()
        print("🎉 All checks passed! gbcms v2 is ready to use.")
        print()
        print("Try running:")
        print("  python -m gbcms.cli --help")
        return 0
    else:
        print()
        print(
            "⚠️  Some checks failed. Please ensure dependencies are installed and the Rust extension is built."
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
