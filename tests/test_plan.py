"""Tests for experiment planner."""

from flex_crispr_probe_designer.plan.experiment import plan_experiment


class TestPlanExperiment:
    def test_small_screen(self) -> None:
        plan = plan_experiment(100)
        assert plan.min_cells == 20_000
        assert plan.wta_hybs == 1
        assert plan.bc_oligo_hybs == 1
        assert plan.gem_wells == 1

    def test_medium_screen(self) -> None:
        plan = plan_experiment(5000)
        # Should round up to the 9600 tier
        assert plan.wta_hybs == 8
        assert plan.gem_wells == 2
        assert plan.bc_oligo_hybs == 48

    def test_large_screen(self) -> None:
        plan = plan_experiment(19000)
        assert plan.wta_hybs == 16
        assert plan.gem_wells == 4

    def test_exact_tier_match(self) -> None:
        plan = plan_experiment(4800)
        assert plan.min_cells == 960_000
        assert plan.wta_hybs == 4
        assert plan.bc_oligo_hybs == 48
        assert plan.gem_wells == 1

    def test_custom_cells_per_guide(self) -> None:
        plan = plan_experiment(100, cells_per_guide=500)
        assert plan.min_cells == 50_000
        assert plan.cells_per_guide == 500

    def test_sequencing_calculations(self) -> None:
        plan = plan_experiment(100)
        assert plan.gex_reads_total == 20_000 * 10_000
        assert plan.crispr_reads_total == 20_000 * 2_000
        assert plan.total_reads == plan.gex_reads_total + plan.crispr_reads_total

    def test_beyond_table_range(self) -> None:
        plan = plan_experiment(40_000)
        # Should extrapolate from largest tier
        assert plan.gem_wells > 4
        assert plan.wta_hybs > 16

    def test_between_tiers_rounds_up(self) -> None:
        # 300 guides is between 100 and 400 tiers
        plan = plan_experiment(300)
        assert plan.wta_hybs == 1
        assert plan.bc_oligo_hybs == 4  # Rounds up to 400 tier
