function UserFun1(obj)
    fprintf(1,'This is a test file.\n After creating this file, it is necessary to restart MATLAB so that MATLAB recognizes this file as a PO method!\n');
    SubUserFun1
end

function SubUserFun1
    fprintf(1,'SubUserFun1 is only accesisble in UserFun1.\n');
end