excel_in = 'K:\My Documents\cricket\scorecards\lge_scorecards.xlsx';

num_games = 22;

%players = [];
teams = cell(11,num_games);
bowlers = [];
batters = [];
fielders = [];
bowlers_stats = struct;
batters_stats = struct;
fielders_stats = struct;

games(num_games).team = [];
games(num_games).bowlers = [];
for game = 1:num_games; 
    
    [d d data] = xlsread(excel_in, game);

    %Workout if we're home or away and who the opposition is
    oppo = data{1,1};
    %home_or_away = 'away';
    if strfind(oppo, 'Hale Barns')
        oppo = data{1,4};
        %home_or_away = 'home';
    end
    %display([home_or_away ' vs ' oppo]);
    
    %Workout whether we batted or bowled
    if strfind(data{3,3}, 'Hale Barns')
        batsman_column = 2;
        bowlers_column = 11;
        catches_column = 12;
    else
        batsman_column = 11;
        bowlers_column = 2;
        catches_column = 3;
    end
    
    %Loop through each batsman adding their stats
    for pp = 1:11
        how_out = data{pp+4, catches_column};
        if isempty(how_out)
            continue;
        elseif strfind(how_out, 'st') == 1
            
            f_name = how_out(4:end);
            f_id = find(strcmp(f_name, fielders));
            
            %If the haven't got one create a new entry
            if isempty(f_id)
                fielders{end+1,1} = f_name; %#ok
                f_id = length(fielders);

                %Initialise stats structure for this bat
                fielders_stats(f_id).name = f_name;
                fielders_stats(f_id).catches = 0;
                fielders_stats(f_id).stumpings = 0;
                fielders_stats(f_id).opposition = [];
                fielders_stats(f_id).how_out = [];
            end
            
            %Increment games played
            fielders_stats(f_id).stumpings = fielders_stats(f_id).stumpings + 1;
            fielders_stats(f_id).opposition{end+1,1} = oppo;
            fielders_stats(f_id).how_out{end+1,1} = 'st';
            
        elseif strfind(how_out, 'ct+b') == 1
            f_name = data{pp+4, catches_column+1}(3:end);
            f_id = find(strcmp(f_name, fielders));
            
            %If the haven't got one create a new entry
            if isempty(f_id)
                fielders{end+1,1} = f_name; %#ok
                f_id = length(fielders);

                %Initialise stats structure for this bat
                fielders_stats(f_id).name = f_name;
                fielders_stats(f_id).catches = 0;
                fielders_stats(f_id).stumpings = 0;
                fielders_stats(f_id).opposition = [];
                fielders_stats(f_id).how_out = [];
            end
            
            %Increment games played
            fielders_stats(f_id).catches = fielders_stats(f_id).catches + 1;
            fielders_stats(f_id).opposition{end+1,1} = oppo;
            fielders_stats(f_id).how_out{end+1,1} = 'ct+b';
            
        elseif strfind(how_out, 'ct') == 1
            f_name = how_out(4:end);
            f_id = find(strcmp(f_name, fielders));
            
            %If the haven't got one create a new entry
            if isempty(f_id)
                fielders{end+1,1} = f_name; %#ok
                f_id = length(fielders);

                %Initialise stats structure for this bat
                fielders_stats(f_id).name = f_name;
                fielders_stats(f_id).catches = 0;
                fielders_stats(f_id).stumpings = 0;
                fielders_stats(f_id).opposition = [];
                fielders_stats(f_id).how_out = [];
            end
            
            %Increment games played
            fielders_stats(f_id).catches = fielders_stats(f_id).catches + 1;
            fielders_stats(f_id).opposition{end+1,1} = oppo;
            fielders_stats(f_id).how_out{end+1,1} = 'ct';
        else
            display([how_out(1:end)]);
        end
    end
end
%%

        %Workout the id for this batsman
        f_id = find(strcmp(data{pp+4, catches_column}, fielders));
        
            %If the haven't got one create a new entry
            if isempty(f_id)
                fielders{end+1,1} = data{pp+4, catches_column}; %#ok
                f_id = length(fielders);

                %Initialise stats structure for this bat
                batters_stats(bat_id).name = batters{bat_id};
                batters_stats(bat_id).played = 0;
                batters_stats(bat_id).innings = 0;
                batters_stats(bat_id).not_outs = 0;
                batters_stats(bat_id).opposition = [];
                batters_stats(bat_id).scores = [];
                batters_stats(bat_id).how_out = [];
            end
        
        %Increment games played
        batters_stats(bat_id).played = batters_stats(bat_id).played + 1;
        
        %Check whether they batted
        score = data{pp+4, batsman_column+5};
        if ~isnan(score)
            
            %Increment innings count and add score and opposition
            batters_stats(bat_id).innings = batters_stats(bat_id).innings + 1;
            batters_stats(bat_id).scores(end+1,1) = score;
            batters_stats(bat_id).opposition{end+1,1} = oppo;
            
            %Check how they were out
            how_out = data{pp+4, batsman_column+1};
            batters_stats(bat_id).how_out{end+1,1} = how_out;
            
            %If not out, increment the not outs score
            if strcmpi(how_out, 'not out')
                batters_stats(bat_id).not_outs = batters_stats(bat_id).not_outs + 1;
            end
        end
        
    end
    
    
    %Workout who bowled for us
    go = true;
    bowler = 1;
    while go
        games(game).bowlers{bowler,1} = data{bowler+23, bowlers_column};
        
        %Workout the id for this batsman
        bowl_id = find(strcmp(data{bowler+23, bowlers_column}, bowlers));
        
        %If the haven't got one create a new entry
        if isempty(bowl_id)
            bowlers{end+1,1} = data{bowler+23, bowlers_column}; %#ok
            bowl_id = length(bowlers);
            
            %Initialise stats structure for this bat
            bowlers_stats(bowl_id).name = bowlers{bowl_id};
            bowlers_stats(bowl_id).opposition = [];
            bowlers_stats(bowl_id).overs = [];
            bowlers_stats(bowl_id).maidens = [];
            bowlers_stats(bowl_id).runs = [];
            bowlers_stats(bowl_id).wickets = [];
        end
        
        %Update runs, wickets etc
        bowlers_stats(bowl_id).opposition{end+1,1} = oppo;
        bowlers_stats(bowl_id).overs(end+1,1) = data{bowler+23, bowlers_column+1};
        bowlers_stats(bowl_id).maidens(end+1,1) = data{bowler+23, bowlers_column+2};
        bowlers_stats(bowl_id).runs(end+1,1) = data{bowler+23, bowlers_column+3};
        bowlers_stats(bowl_id).wickets(end+1,1) = data{bowler+23, bowlers_column+4};
        
        %See if anyone else bowled in this match
        bowler = bowler + 1;
        go = ~isnan(data{bowler+23, bowlers_column});
    end
    %
end

num_bowlers = length(bowlers);
num_batsman = length(batters);

if ~isempty(excel_out)
    
    %Make containers to add data to excel file
    bowl_overs = zeros(num_bowlers, 1);
    bowl_mdns = zeros(num_bowlers, 1);
    bowl_runs = zeros(num_bowlers, 1);
    bowl_wkts = zeros(num_bowlers, 1);
    bowl_averages = zeros(num_bowlers, 1);
    bowl_srs = zeros(num_bowlers, 1);
    bowl_ers = zeros(num_bowlers, 1);
    bowl_bfs = cell(num_bowlers, 1);

    bat_games = zeros(num_batsman,1);
    bat_inns = zeros(num_batsman,1);
    bat_nos = zeros(num_batsman,1);
    bat_runs = zeros(num_batsman,1);
    bat_averages = zeros(num_batsman,1);
    bat_hs = zeros(num_batsman,1);
    
    %Make excel file headers
    xlswrite(excel_out, {'Name'}, 'Batting', 'A1');
    xlswrite(excel_out, {'Games'}, 'Batting', 'B1');
    xlswrite(excel_out, {'Inns'}, 'Batting', 'C1');
    xlswrite(excel_out, {'Not outs'}, 'Batting', 'D1');
    xlswrite(excel_out, {'Runs'}, 'Batting', 'E1');
    xlswrite(excel_out, {'Average'}, 'Batting', 'F1');
    xlswrite(excel_out, {'HS'}, 'Batting', 'G1');
    
    xlswrite(excel_out, {'Name'}, 'Bowling', 'A1');
    xlswrite(excel_out, {'Overs'}, 'Bowling', 'B1');
    xlswrite(excel_out, {'Mdns'}, 'Bowling', 'C1');
    xlswrite(excel_out, {'Runs'}, 'Bowling', 'D1');
    xlswrite(excel_out, {'Wkts'}, 'Bowling', 'E1');
    xlswrite(excel_out, {'Average'}, 'Bowling', 'F1');
    xlswrite(excel_out, {'Strike Rate'}, 'Bowling', 'G1');
    xlswrite(excel_out, {'Econ Rate'}, 'Bowling', 'H1');
    xlswrite(excel_out, {'Best Figs'}, 'Bowling', 'I1');
end


%Get cumulative stats for the batsmen
for bat_id = 1:num_batsman
    batters_stats(bat_id).runs_scored = sum(batters_stats(bat_id).scores);
    batters_stats(bat_id).average = batters_stats(bat_id).runs_scored / ...
        (batters_stats(bat_id).innings - batters_stats(bat_id).not_outs);
    batters_stats(bat_id).highest_score = max(batters_stats(bat_id).scores);
    
    if ~isempty(excel_out)
        %Add stats to the containers for the excel spreadsheet general
        %stats
        bat_games(bat_id) = batters_stats(bat_id).played;
        bat_inns(bat_id) = batters_stats(bat_id).innings;
        bat_nos(bat_id) = batters_stats(bat_id).not_outs;
        bat_runs(bat_id) = batters_stats(bat_id).runs_scored;
        if ~isempty(batters_stats(bat_id).highest_score)
            bat_averages(bat_id) = batters_stats(bat_id).average;
            bat_hs(bat_id) = batters_stats(bat_id).highest_score; 
        end
        
        %Now write out the personal stats
        end_cell = num2str(batters_stats(bat_id).innings + 5);
        xlswrite(excel_out, {'Batting'}, batters{bat_id}, 'A1');
        xlswrite(excel_out, batters(bat_id), batters{bat_id}, 'B1');
        
        xlswrite(excel_out, {'Games'}, batters{bat_id}, 'A2');
        xlswrite(excel_out, {'Inns'}, batters{bat_id}, 'B2');
        xlswrite(excel_out, {'Not outs'}, batters{bat_id}, 'C2');
        xlswrite(excel_out, {'Runs'}, batters{bat_id}, 'D2');
        xlswrite(excel_out, {'Average'}, batters{bat_id}, 'E2');
        xlswrite(excel_out, {'HS'}, batters{bat_id}, 'F2');
        
        xlswrite(excel_out, batters_stats(bat_id).played, batters{bat_id}, 'A3');
        xlswrite(excel_out, batters_stats(bat_id).innings, batters{bat_id}, 'B3');
        xlswrite(excel_out, batters_stats(bat_id).not_outs, batters{bat_id}, 'C3');
        xlswrite(excel_out, batters_stats(bat_id).runs_scored, batters{bat_id}, 'D3');
        if ~isempty(batters_stats(bat_id).highest_score)
            xlswrite(excel_out, batters_stats(bat_id).average, batters{bat_id}, 'E3');
            xlswrite(excel_out, batters_stats(bat_id).highest_score, batters{bat_id}, 'F3');
            
            xlswrite(excel_out, {'Opposition'}, batters{bat_id}, 'A5');
            xlswrite(excel_out,  batters_stats(bat_id).opposition, batters{bat_id}, ['A6:A' end_cell]);

            xlswrite(excel_out, {'Runs'}, batters{bat_id}, 'B5');
            xlswrite(excel_out,  batters_stats(bat_id).scores, batters{bat_id}, ['B6:B' end_cell]);

            xlswrite(excel_out, {'How Out'}, batters{bat_id}, 'C5');
            xlswrite(excel_out,  batters_stats(bat_id).how_out, batters{bat_id}, ['C6:C' end_cell]);
        end   
    end
   
   
end

