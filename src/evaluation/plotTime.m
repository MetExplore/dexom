% Script to export the time of each method used for the Yeast 6 evaluation
% Note that this scripts takes a zip file with the results and generates
% the CSV. The evaluation module contains the required zip file in
% evaluation/yeast/essential_genes/yeast6-evaluation/matlab/dexom-yeast6-eval-matlab.zip
% This file is a multipart-zip file that should be converted to a single
% zip file as Matlab cannot extract these types of zip files.

% Exported files are located in the dexom-yeast6-eval-export-csv.zip file

tmp = tempname();
mkdir(tmp);
% Name of the single zip file with all the results
unzip('yeast6-evaluation-nogit.zip', tmp);
files = dir([tmp filesep '07-06-2020/*-th-0.25-0.75-*.mat']);
figure;
hold on;
times = [];
solutions = [];
labels = [];
for i = 1:numel(files)
    f = [files(i).folder filesep files(i).name];
    fprintf('Ploting %s...\n', f);
    % Load the result file and aggregate the essential genes
    s = load(f, 'store');
    time = s.store.result.elapsedTime;
    sols = cumsum(s.store.result.accepted);
    plot(time/60, sols, 'DisplayName', s.store.info.method);
    times = [times; time'/60];
    solutions = [solutions; sols'];
    labels = [labels; repelem(string(s.store.info.method), numel(sols))'];
end
hold off;
legend;

t = table(times, solutions, labels);
writetable(t, 'time_th-0.25-0.75.csv');