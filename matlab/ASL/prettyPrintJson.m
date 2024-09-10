function prettyJson = prettyPrintJson(jsonString)
    indent = '    ';  % Define the indentation (4 spaces)
    level = 0;  % Keep track of the indentation level
    prettyJson = '';  % Initialize the pretty-printed JSON string
    
    i = 1;
    insideString = false;  % Track whether we are inside a string

    while i <= length(jsonString)
        char = jsonString(i);
        
        switch char
            case '"'
                % Toggle insideString flag when we encounter a quote
                insideString = ~insideString;
                prettyJson = append(prettyJson, char);
                
            case '{'
                if ~insideString
                    level = level + 1;
                    prettyJson = append(prettyJson, char, newline, repmat(indent, 1, level));
                else
                    prettyJson = append(prettyJson, char);
                end
                
            case '['
                if ~insideString
                    level = level + 1;
                    prettyJson = append(prettyJson, char, newline, repmat(indent, 1, level));
                else
                    prettyJson = append(prettyJson, char);
                end
                
            case '}'
                if ~insideString
                    level = level - 1;
                    prettyJson = append(prettyJson, newline, repmat(indent, 1, level), char);
                else
                    prettyJson = append(prettyJson, char);
                end
                
            case ']'
                if ~insideString
                    level = level - 1;
                    prettyJson = append(prettyJson, newline, repmat(indent, 1, level), char);
                else
                    prettyJson = append(prettyJson, char);
                end
                
            case ','
                if ~insideString
                    prettyJson = append(prettyJson, char, newline, repmat(indent, 1, level));
                else
                    prettyJson = append(prettyJson, char);
                end
                
            case ':'
                if ~insideString
                    prettyJson = append(prettyJson, char, ' ');
                else
                    prettyJson = append(prettyJson, char);
                end
                
            otherwise
                % Escape special characters like newline (\n) if inside a string
                if insideString && char == newline
                    prettyJson = append(prettyJson, '\\n');
                else
                    prettyJson = append(prettyJson, char);
                end
        end
        
        i = i + 1;
    end
end
