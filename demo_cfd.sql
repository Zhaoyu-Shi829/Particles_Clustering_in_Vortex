CREATE TABLE cfd_simulations (
    simulation_id INT PRIMARY KEY,
    date_run TIMESTAMP,
    model_name VARCHAR(100),
    reynolds_number FLOAT,
    mach_number FLOAT,
    inlet_velocity FLOAT,
    outlet_pressure FLOAT,
    num_iterations INT,
    convergence_status VARCHAR(20),
    max_pressure FLOAT,
    avg_velocity FLOAT,
    drag_coefficient FLOAT,
    lift_coefficient FLOAT,
    raw_data_path VARCHAR(255)
);

-- Example insertion
INSERT INTO cfd_simulations 
(simulation_id, date_run, model_name, reynolds_number, mach_number, inlet_velocity, 
 outlet_pressure, num_iterations, convergence_status, max_pressure, avg_velocity, 
 drag_coefficient, lift_coefficient, raw_data_path)
VALUES 
(1, '2024-10-14 10:00:00', 'airfoil_model_v1', 1000000, 0.3, 100, 
 101325, 1000, 'converged', 150000, 80, 0.05, 0.8, 
 '/data/cfd/sim_1_results.h5');

-- Example query
SELECT model_name, reynolds_number, drag_coefficient
FROM cfd_simulations
WHERE mach_number < 0.5 AND convergence_status = 'converged'
ORDER BY drag_coefficient ASC
LIMIT 10;
