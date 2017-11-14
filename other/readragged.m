function a = readragged(fn)
ca = {};
fid = fopen(fn, 'r');
while 1
    temp = fgetl(fid);
    if (temp==-1)
        break;
    end
    A = sscanf(temp, '%g');
    ca{end+1} = A;
end
fclose(fid);

fixed = true;
fl = length(ca{1});
for i=1:length(ca)
    if (length(ca{i}) ~= fl)
        fixed = false;
        break;
    end
end
a = ca;
if fixed
    cl = class(ca{1}(1));
    a = zeros(length(ca), fl, cl);
    for i=1:length(ca)
        a(i,:) = ca{i};
    end
end
end