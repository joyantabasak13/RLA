from splink.duckdb.linker import DuckDBLinker
import splink.duckdb.comparison_library as cl
import splink.duckdb.comparison_template_library as ctl
import splink.duckdb.blocking_rule_library as brl
from splink.datasets import splink_datasets

df = splink_datasets.fake_1000

settings = {
    "link_type": "dedupe_only",
    "blocking_rules_to_generate_predictions": [
        brl.exact_match_rule("first_name"),
        brl.exact_match_rule("surname"),
    ],
    "comparisons": [
        ctl.name_comparison("first_name"),
        ctl.name_comparison("surname"),
        ctl.date_comparison("dob", cast_strings_to_date=True),
        cl.exact_match("city", term_frequency_adjustments=True),
        ctl.email_comparison("email", include_username_fuzzy_level=False),
    ],
}

linker = DuckDBLinker(df, settings)
linker.estimate_u_using_random_sampling(max_pairs=1e6)

blocking_rule_for_training = brl.and_(
                            brl.exact_match_rule("first_name"),
                            brl.exact_match_rule("surname")
                            )

linker.estimate_parameters_using_expectation_maximisation(blocking_rule_for_training)

blocking_rule_for_training = brl.exact_match_rule("dob")
linker.estimate_parameters_using_expectation_maximisation(blocking_rule_for_training)


pairwise_predictions = linker.predict()

clusters = linker.cluster_pairwise_predictions_at_threshold(pairwise_predictions, 0.95)
clusters.as_pandas_dataframe(limit=5)
