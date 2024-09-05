function prettyJson = prettyPrintJson(jsonString)
    indent = '    ';  % Define the indentation (4 spaces)
    level = 0;  % Keep track of the indentation level
    prettyJson = '';  % Initialize the pretty-printed JSON string
    
    i = 1;
    while i <= length(jsonString)
        char = jsonString(i);
        
        switch char
            case '{'
                level = level + 1;
                prettyJson = append(prettyJson, char, newline, repmat(indent, 1, level));
                
            case '['
                level = level + 1;
                prettyJson = append(prettyJson, char, newline, repmat(indent, 1, level));
                
            case '}'
                level = level - 1;
                prettyJson = append(prettyJson, newline, repmat(indent, 1, level), char);
                
            case ']'
                level = level - 1;
                prettyJson = append(prettyJson, newline, repmat(indent, 1, level), char);
                
            case ','
                prettyJson = append(prettyJson, char, newline, repmat(indent, 1, level));
                
            case ':'
                prettyJson = append(prettyJson, char, ' ');
                
            otherwise
                prettyJson = append(prettyJson, char);
        end
        
        i = i + 1;
    end
end
