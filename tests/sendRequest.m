% Send HTTP request, get response and additional information in retInfos.
function [response, retInfos, history] = sendRequest(uri, request)

% uri: matlab.net.URI
% request: matlab.net.http.RequestMessage
% response: matlab.net.http.ResponseMessage

% matlab.net.http.HTTPOptions persists across requests to reuse  previous
% Credentials in it for subsequent authentications
persistent options 

% infos is a containers.Map object where: 
%    key is uri.Host; 
%    value is "info" struct containing:
%        cookies: vector of matlab.net.http.Cookie or empty
%        uri: target matlab.net.URI if redirect, or empty
persistent infos

if isempty(options)
    options = matlab.net.http.HTTPOptions('ConnectTimeout',20);
end

if isempty(infos)
    infos = containers.Map;
end
host = string(uri.Host); % get Host from URI
try
    % get info struct for host in map
    info = infos(char(host));
    if ~isempty(info.uri)
        % If it has a uri field, it means a redirect previously
        % took place, so replace requested URI with redirect URI.
        uri = info.uri;
    end
    if ~isempty(info.cookies)
        % If it has cookies, it means we previously received cookies from this host.
        % Add Cookie header field containing all of them.
        request = request.addFields(matlab.net.http.field.CookieField(info.cookies));
    end
catch
    % no previous redirect or cookies for this host
    info = [];
end

% Send request and get response and history of transaction.
[response, ~, history] = request.send(uri, options);
if response.StatusCode ~= matlab.net.http.StatusCode.OK
    return
end

% Get the Set-Cookie header fields from response message in
% each history record and save them in the map.
arrayfun(@addCookies, history)

% If the last URI in the history is different from the URI sent in the original 
% request, then this was a redirect. Save the new target URI in the host info struct.
targetURI = history(end).URI;
if ~isequal(targetURI, uri)
    if isempty(info)
        % no previous info for this host in map, create new one
        infos(char(host)) = struct('cookies',[],'uri',targetURI);
    else
        % change URI in info for this host and put it back in map
        info.uri = targetURI;
        infos(char(host)) = info;
    end
end
retInfos = infos;

    function addCookies(record)
        % Add cookies in Response message in history record
        % to the map entry for the host to which the request was directed.
        %
        ahost = record.URI.Host; % the host the request was sent to
        cookieFields = record.Response.getFields('Set-Cookie');
        if isempty(cookieFields)
            return
        end
        cookieData = cookieFields.convert(); % get array of Set-Cookie structs
        cookies = [cookieData.Cookie]; % get array of Cookies from all structs
        try
            % If info for this host was already in the map, add its cookies to it.
            ainfo = infos(ahost);
            ainfo.cookies = [ainfo.cookies cookies];
            infos(char(ahost)) = ainfo;
        catch
            % Not yet in map, so add new info struct.
            infos(char(ahost)) = struct('cookies',cookies,'uri',[]);
        end
    end
end