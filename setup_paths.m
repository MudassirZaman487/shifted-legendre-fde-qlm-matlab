% Add research-code directories from the repository root.
repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'src', 'operators'));
addpath(fullfile(repoRoot, 'src', 'linear_solvers'));
addpath(fullfile(repoRoot, 'src', 'qlm'));
addpath(fullfile(repoRoot, 'examples', 'thesis'));
addpath(fullfile(repoRoot, 'examples', 'matrix_demos'));
addpath(fullfile(repoRoot, 'tests'));
disp('Shifted Legendre thesis MATLAB paths added.');
