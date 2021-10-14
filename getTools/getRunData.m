function s = getRunData(folder)
    pingName = dir([folder 'data*']);
    s = load([folder pingName.name]);
end