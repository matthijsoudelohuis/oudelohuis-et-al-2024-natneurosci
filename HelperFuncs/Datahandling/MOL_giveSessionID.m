function session_ID = MOL_giveSessionID(mousename,rec_datetime,expnumber)

%Get the mouse part:
if strfind(mousename,'.')
    splitmouse_ID = strsplit(mousename,'.');
    if strfind(splitmouse_ID{1},'P')
        splitmouse_ID{1} = splitmouse_ID{1}(2);
    end
    mouse = sprintf('%01d%02d',str2double(splitmouse_ID{1}),str2double(splitmouse_ID{2}));
else
    mouse = mousename;
end

%Get the datetime part:
if strfind(rec_datetime,'_')
    splitrec_datetime = strsplit(rec_datetime,'_');
    splitrec_datetime = strsplit(splitrec_datetime{1},'-');
    session = sprintf('%04d%02d%02d',str2double(splitrec_datetime{1}),str2double(splitrec_datetime{2}),str2double(splitrec_datetime{3}));
else
    session = rec_datetime;
    if numel(session)==12
        session = [session '00'];
    end
end

%Get the experiment number:
expnumber = num2str(expnumber);
switch numel(expnumber)
    case 1
        expnumber = ['0' expnumber];
    case 2
    case {3,4,5,6}
        expnumber = expnumber(end-1:end);
end

%Compile to unique session_ID
session_ID = [mouse session expnumber];

end