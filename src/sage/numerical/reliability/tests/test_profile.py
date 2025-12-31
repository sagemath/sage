from sage.numerical.reliability.profile import ReliabilityProfile


def test_default_profile_state():
    profile = ReliabilityProfile()
    summary = profile.summary()

    assert summary["trust_level"] == "UNKNOWN"
    assert summary["failure_signals"] == []
    assert summary["assumptions"] == []


def test_profile_updates():
    profile = ReliabilityProfile()

    profile.set_trust_level("HIGH")
    profile.add_failure_signal("LOSS_OF_SIGNIFICANCE")
    profile.add_assumption("x > 0")

    summary = profile.summary()

    assert summary["trust_level"] == "HIGH"
    assert "LOSS_OF_SIGNIFICANCE" in summary["failure_signals"]
    assert "x > 0" in summary["assumptions"]
