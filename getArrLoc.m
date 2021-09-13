function loc = getArrLoc(hID)

switch lower(hID)
    case 'av'
        loc = [49.0907 64.5372];
    case 'cd'
        loc = [49.1933  64.8315];
    case 'mb'
        loc = [48.6039 64.1863];
    case 'pc'
        loc= [48.5311 64.1983];
    otherwise
        error(['Can''t location. AV = Anse-a-Valleau, CD = Cloridorme, MB = Malbay, PC = Perce.' ])
end
        