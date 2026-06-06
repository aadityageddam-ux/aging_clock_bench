"""AgingClockBench: Benchmark biological aging clocks on your data."""

from agingclockbench.clocks.phenoage import PhenoAge
from agingclockbench.clocks.kdm import KDM
from agingclockbench.clocks.dunedinpace import DunedinPACEProxy
from agingclockbench.benchmarks.suite import BenchmarkSuite

__version__ = "0.1.0"
__all__ = ["PhenoAge", "KDM", "DunedinPACEProxy", "BenchmarkSuite"]
