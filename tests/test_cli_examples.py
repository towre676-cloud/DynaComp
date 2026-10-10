import json
import unittest
from pathlib import Path
from dynacomp.cli import compile_spec, replay

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


class ExampleSpecTests(unittest.TestCase):
    def test_every_example_compiles_round_trips_and_replays(self):
        specs = sorted(EXAMPLES.glob("*.json"))
        self.assertGreaterEqual(len(specs), 11)
        for path in specs:
            with self.subTest(example=path.name):
                packet = json.loads(json.dumps(compile_spec(json.loads(path.read_text()))))
                self.assertTrue(replay(packet))

    def test_unknown_schema_rejected(self):
        with self.assertRaises(ValueError):
            replay({"schema": "dynacomp.unknown/1"})
