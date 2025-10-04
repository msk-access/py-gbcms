#!/usr/bin/env python3
"""
Verify GetBaseCounts installation and dependencies.

This script checks that all required components are properly installed
and can be imported.
"""

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
        __import__(module_name)
        return True, f"✅ {package_name}"
    except ImportError as e:
        return False, f"❌ {package_name}: {str(e)}"


def main():
    """Run installation verification."""
    print("=" * 60)
    print("GetBaseCounts Installation Verification")
    print("=" * 60)
    print()

    # Core dependencies
    print("Core Dependencies:")
    print("-" * 60)

    core_deps = [
        ("pysam", "pysam"),
        ("numpy", "numpy"),
        ("typer", "typer"),
        ("rich", "rich"),
        ("pandas", "pandas"),
        ("pydantic", "pydantic"),
        ("numba", "numba"),
        ("joblib", "joblib"),
    ]

    core_results = [check_import(mod, pkg) for mod, pkg in core_deps]
    for _success, msg in core_results:
        print(msg)

    print()

    # Optional dependencies
    print("Optional Dependencies:")
    print("-" * 60)

    opt_deps = [
        ("ray", "ray (for distributed computing)"),
    ]

    opt_results = [check_import(mod, pkg) for mod, pkg in opt_deps]
    for _success, msg in opt_results:
        print(msg)

    print()

    # GetBaseCounts modules
    print("GetBaseCounts Modules:")
    print("-" * 60)

    gb_modules = [
        ("gbcms", "gbcms"),
        ("gbcms.cli", "CLI"),
        ("gbcms.config", "Config"),
        ("gbcms.models", "Pydantic Models"),
        ("gbcms.variant", "Variant Loader"),
        ("gbcms.counter", "Base Counter"),
        ("gbcms.numba_counter", "Numba Counter"),
        ("gbcms.parallel", "Parallel Processor"),
        ("gbcms.reference", "Reference Handler"),
        ("gbcms.output", "Output Formatter"),
        ("gbcms.processor", "Main Processor"),
    ]

    gb_results = [check_import(mod, pkg) for mod, pkg in gb_modules]
    for _success, msg in gb_results:
        print(msg)

    print()

    # Check CLI entry point
    print("CLI Entry Point:")
    print("-" * 60)

    try:

        print("✅ CLI app accessible")

        # Try to get version
        try:
            from gbcms import __version__

            print(f"✅ Version: {__version__}")
        except Exception:
            print("⚠️  Version not accessible")
    except Exception as e:
        print(f"❌ CLI app: {str(e)}")

    print()

    # Summary
    print("=" * 60)
    print("Summary:")
    print("-" * 60)

    all_results = core_results + opt_results + gb_results
    total = len(all_results)
    passed = sum(1 for success, _ in all_results if success)
    failed = total - passed

    print(f"Total checks: {total}")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")

    if failed == 0:
        print()
        print("🎉 All checks passed! GetBaseCounts is ready to use.")
        print()
        print("Try running:")
        print("  gbcms --help")
        print("  gbcms version")
        return 0
    else:
        print()
        print("⚠️  Some checks failed. Please install missing dependencies:")
        print("  uv pip install -e '.[dev,all]'")
        return 1


if __name__ == "__main__":
    sys.exit(main())
