PRAGMA foreign_keys = ON;

CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    condition TEXT NOT NULL,
    age INTEGER NOT NULL CHECK (age >= 0),
    sex TEXT NOT NULL CHECK (sex IN ('M', 'F')),
    treatment TEXT NOT NULL,
    response TEXT CHECK (response IN ('yes', 'no') OR response IS NULL)
);

CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY,
    subject_id TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL
);

CREATE TABLE cell_populations (
    population TEXT PRIMARY KEY
);

CREATE TABLE cell_counts (
    sample_id TEXT NOT NULL REFERENCES samples(sample_id),
    population TEXT NOT NULL REFERENCES cell_populations(population),
    count INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population)
);

CREATE INDEX idx_subject_cohort ON subjects(condition, treatment, response);
CREATE INDEX idx_samples_time_type ON samples(time_from_treatment_start, sample_type);
CREATE INDEX idx_samples_subject ON samples(subject_id);
CREATE INDEX idx_cells_population ON cell_counts(population);

CREATE VIEW sample_frequencies AS
WITH totals AS (
    SELECT sample_id, SUM(count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
)
SELECT cc.sample_id AS sample,
       t.total_count AS total_count,
       cc.population AS population,
       cc.count AS count,
       (100.0 * cc.count / NULLIF(t.total_count, 0)) AS percentage
FROM cell_counts cc
JOIN totals t ON cc.sample_id = t.sample_id;

CREATE VIEW frequency_details AS
SELECT sf.sample, sf.total_count, sf.population, sf.count, sf.percentage,
       s.subject_id, s.sample_type, s.time_from_treatment_start,
       su.project_id AS project, su.condition, su.treatment, su.response, su.sex
FROM sample_frequencies sf
JOIN samples s ON s.sample_id = sf.sample
JOIN subjects su ON su.subject_id = s.subject_id;
