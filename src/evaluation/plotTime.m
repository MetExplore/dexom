tmp = tempname();
mkdir(tmp);
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