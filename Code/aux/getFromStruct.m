function myval = getFromStruct(mystruct,myfield,defaultval)
    % get the value of the specified field from a struct if the field exists,
    % otherwise set to default value
    %
    % INPUTS:
    % - mystruct: a struct variable
    % - myfield: a string (supposed to be the field name)
    % - defaultval: default value to set when field is nonexisting
    % RETURNS:
    % - myval: the value (either from the struct or the specified default)

    if(isfield(mystruct,myfield))
        myval = mystruct.(myfield);
    else
        myval = defaultval; % set default value
    end

end
