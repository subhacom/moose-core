"""Structured account of what a load did and did not cover.

The whole point of the new reader is to never silently mis-simulate: whatever
we cannot map faithfully must be visible. A ``LoadReport`` is returned (and
attached to the handler) so callers can inspect coverage programmatically and
so tests can track it as a regression metric.
"""

from dataclasses import dataclass, field


@dataclass
class LoadReport:
    filepath: str = ''
    loadpath: str = ''

    # counts of what mapped, and how
    reactions_native: int = 0      # Reac / MMenz
    reactions_function: int = 0    # arbitrary rate law via Function
    rate_rules: int = 0
    assignment_rules: int = 0

    normalized: list = field(default_factory=list)

    # things we could not represent (each entry is a human-readable reason)
    unsupported: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    def unsupported_add(self, msg):
        self.unsupported.append(msg)

    @property
    def fully_supported(self):
        return not self.unsupported

    def summary(self):
        lines = [
            'SBML load report for %s -> %s' % (self.filepath, self.loadpath),
            '  reactions: %d native, %d via Function'
            % (self.reactions_native, self.reactions_function),
            '  rules: %d rate, %d assignment'
            % (self.rate_rules, self.assignment_rules),
        ]
        if self.unsupported:
            lines.append('  UNSUPPORTED (%d):' % len(self.unsupported))
            lines += ['    - ' + u for u in self.unsupported]
        if self.warnings:
            lines.append('  warnings (%d):' % len(self.warnings))
            lines += ['    - ' + w for w in self.warnings]
        return '\n'.join(lines)

    def __str__(self):
        return self.summary()
