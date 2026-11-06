function ids = rsa_image_ids(values)
%RSA_IMAGE_IDS Normalize image1.png through image400.png stimulus labels.
ids = stems(values);
assert(all(ismember(ids,"image"+string((1:400)'))), ...
    'Stimulus labels must use image1.png through image400.png.');
end

function result = stems(values)
values = string(values(:)); result = strings(size(values));
for i = 1:numel(values)
    [~, name, ~] = fileparts(strrep(char(values(i)), '\', '/'));
    result(i) = lower(strtrim(string(name)));
end
end
