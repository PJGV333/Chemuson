"""Pruebas de política de chequeo de updates."""

from datetime import datetime, timedelta, timezone


from chemuson.update.policy import (
    can_offer_update,
    is_downgrade_suspected,
    mark_checked,
    should_check_now,
)
from chemuson.update.types import UpdateSettings


def test_should_check_now_respects_interval() -> None:
    now = datetime(2026, 2, 26, 12, 0, tzinfo=timezone.utc)
    settings = UpdateSettings(enabled=True, check_interval_hours=24)
    assert should_check_now(settings, now=now) is True
    mark_checked(settings, now=now)
    assert should_check_now(settings, now=now + timedelta(hours=1)) is False
    assert should_check_now(settings, now=now + timedelta(hours=24, minutes=1)) is True


def test_disabled_policy_never_checks() -> None:
    settings = UpdateSettings(enabled=False, check_interval_hours=1)
    assert should_check_now(settings) is False


def test_policy_channel_rules_stable_vs_beta() -> None:
    assert can_offer_update(
        channel="stable",
        current_version="1.0.0",
        candidate_version="1.1.0",
    ) is True
    assert can_offer_update(
        channel="stable",
        current_version="1.0.0",
        candidate_version="1.1.0-beta.1",
    ) is False
    assert can_offer_update(
        channel="beta",
        current_version="1.0.0",
        candidate_version="1.1.0-beta.1",
    ) is True


def test_policy_invalid_last_check_falls_back_to_check_now() -> None:
    settings = UpdateSettings(enabled=True, check_interval_hours=24, last_check_iso="not-iso")
    assert should_check_now(settings) is True


def test_policy_detects_downgrade_replay() -> None:
    assert is_downgrade_suspected("1.4.0", "1.5.0") is True
    assert is_downgrade_suspected("1.5.0", "1.5.0") is False
    assert is_downgrade_suspected("1.6.0", "1.5.0") is False
