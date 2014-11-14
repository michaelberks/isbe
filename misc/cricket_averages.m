function [games bowlers bowlers_stats batters batters_stats data] = cricket_averages(excel_in, excel_out, num_games)

%players = [];
teams = cell(11,num_games);
bowlers = [];
batters = [];
bowlers_stats = struct;
batters_stats = struct;

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
    else
        batsman_column = 11;
        bowlers_column = 2;
    end
    
    %Workout who played for us and batting order
    games(game).team = data(5:15, batsman_column);
    if length(unique(games(game).team)) ~= 11
        display(games(game).team);
        error('Duplicate name in team');
    end
    
    %Loop through each batsman adding their stats
    for pp = 1:11
        teams{pp, game} = data{pp+4, batsman_column};
        
        %Workout the id for this batsman
        bat_id = find(strcmp(data{pp+4, batsman_column}, batters));
        
        %If the haven't got one create a new entry
        if isempty(bat_id)
            batters{end+1,1} = data{pp+4, batsman_column}; %#ok
            bat_id = length(batters);
            
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

%Get cumulative stats for the bowlers
for bowl_id = 1:num_bowlers
    
    %For overs sum we need to account for incomplete overs
    complete_overs = floor(bowlers_stats(bowl_id).overs);
    extra_balls = sum(10*(bowlers_stats(bowl_id).overs - complete_overs));
    total_balls = 6*sum(complete_overs) + extra_balls;
    
    extra_overs = floor(extra_balls / 6);
    extra_balls = rem(extra_balls, 6) / 10;
    
    bowlers_stats(bowl_id).overs_total = sum(complete_overs) + extra_overs + extra_balls;
    
    %For runs, maidens, wickets just a simple sum
    bowlers_stats(bowl_id).maidens_total = sum(bowlers_stats(bowl_id).maidens);
    bowlers_stats(bowl_id).runs_total = sum(bowlers_stats(bowl_id).runs);
    bowlers_stats(bowl_id).wickets_total = sum(bowlers_stats(bowl_id).wickets);
    
    %Compute average and strike rate
    bowlers_stats(bowl_id).average = bowlers_stats(bowl_id).runs_total / ...
        bowlers_stats(bowl_id).wickets_total;
    bowlers_stats(bowl_id).strike_rate = total_balls / ...
        bowlers_stats(bowl_id).wickets_total;
    bowlers_stats(bowl_id).economy_rate = (bowlers_stats(bowl_id).runs_total / total_balls)*6;
    
    %Workout best figures
    most_wickets = max(bowlers_stats(bowl_id).wickets);
    least_runs = min(bowlers_stats(bowl_id).runs(bowlers_stats(bowl_id).wickets == most_wickets));
    bowlers_stats(bowl_id).best_figures = [num2str(most_wickets) '/' num2str(least_runs)];
    
    if ~isempty(excel_out)
        %Add to the containers for the excel spreadsheet
        bowl_overs(bowl_id) = bowlers_stats(bowl_id).overs_total;
        bowl_mdns(bowl_id) = bowlers_stats(bowl_id).maidens_total;
        bowl_runs(bowl_id) = bowlers_stats(bowl_id).runs_total;
        bowl_wkts(bowl_id) = bowlers_stats(bowl_id).wickets_total;
        bowl_averages(bowl_id) = bowlers_stats(bowl_id).average;
        bowl_srs(bowl_id) = bowlers_stats(bowl_id).strike_rate;
        bowl_ers(bowl_id) = bowlers_stats(bowl_id).economy_rate;
        bowl_bfs{bowl_id} = bowlers_stats(bowl_id).best_figures;
        
        xlswrite(excel_out, {'Bowling'}, bowlers{bowl_id}, 'H1');
        xlswrite(excel_out, {'Overs'}, bowlers{bowl_id}, 'H2');
        xlswrite(excel_out, {'Mdns'}, bowlers{bowl_id}, 'I2');
        xlswrite(excel_out, {'Runs'}, bowlers{bowl_id}, 'J2');
        xlswrite(excel_out, {'Wkts'}, bowlers{bowl_id}, 'K2');
        xlswrite(excel_out, {'Average'}, bowlers{bowl_id}, 'L2');
        xlswrite(excel_out, {'Strike Rate'}, bowlers{bowl_id}, 'M2');
        xlswrite(excel_out, {'Econ Rate'}, bowlers{bowl_id}, 'N2');
        xlswrite(excel_out, {'Best Figs'}, bowlers{bowl_id}, 'O2');
        
        xlswrite(excel_out, bowlers_stats(bowl_id).overs_total, bowlers{bowl_id}, 'H3');
        xlswrite(excel_out, bowlers_stats(bowl_id).maidens_total, bowlers{bowl_id}, 'I3');
        xlswrite(excel_out, bowlers_stats(bowl_id).runs_total, bowlers{bowl_id}, 'J3');
        xlswrite(excel_out, bowlers_stats(bowl_id).wickets_total, bowlers{bowl_id}, 'K3');
        xlswrite(excel_out, bowlers_stats(bowl_id).average, bowlers{bowl_id}, 'L3');
        xlswrite(excel_out, bowlers_stats(bowl_id).strike_rate, bowlers{bowl_id}, 'M3');
        xlswrite(excel_out, bowlers_stats(bowl_id).economy_rate, bowlers{bowl_id}, 'N3');
        xlswrite(excel_out, {bowlers_stats(bowl_id).best_figures}, bowlers{bowl_id}, 'O3');
        
        xlswrite(excel_out, {'Opposition'}, bowlers{bowl_id}, 'H5');
        xlswrite(excel_out, {'Overs'}, bowlers{bowl_id}, 'I5');
        xlswrite(excel_out, {'Mdns'}, bowlers{bowl_id}, 'J5');
        xlswrite(excel_out, {'Runs'}, bowlers{bowl_id}, 'K5');
        xlswrite(excel_out, {'Wkts'}, bowlers{bowl_id}, 'L5');
        
        end_cell = num2str(length(bowlers_stats(bowl_id).overs) + 5); 
        xlswrite(excel_out, bowlers_stats(bowl_id).opposition, bowlers{bowl_id}, ['H6:H' end_cell]);
        xlswrite(excel_out, bowlers_stats(bowl_id).overs, bowlers{bowl_id}, ['I6:I' end_cell]);
        xlswrite(excel_out, bowlers_stats(bowl_id).maidens, bowlers{bowl_id}, ['J6:J' end_cell]);
        xlswrite(excel_out, bowlers_stats(bowl_id).runs, bowlers{bowl_id}, ['K6:K' end_cell]);
        xlswrite(excel_out, bowlers_stats(bowl_id).wickets, bowlers{bowl_id}, ['L6:L' end_cell]);

        
    end
    
end

if ~isempty(excel_out)
    %Write out the stats to the excel file
    end_cell = num2str(num_batsman+1);
    xlswrite(excel_out, batters, 'Batting', ['A2:A' end_cell]);
    xlswrite(excel_out, bat_games, 'Batting', ['B2:B' end_cell]);
    xlswrite(excel_out, bat_inns, 'Batting', ['C2:C' end_cell]);
    xlswrite(excel_out, bat_nos, 'Batting', ['D2:D' end_cell]);
    xlswrite(excel_out, bat_runs, 'Batting', ['E2:E' end_cell]);
    xlswrite(excel_out, bat_averages, 'Batting', ['F2:F' end_cell]);
    xlswrite(excel_out, bat_hs, 'Batting', ['G2:G' end_cell]);

    end_cell = num2str(num_bowlers+1);
    xlswrite(excel_out, bowlers, 'Bowling', ['A2:A' end_cell]);
    xlswrite(excel_out, bowl_overs, 'Bowling', ['B2:B' end_cell]);
    xlswrite(excel_out, bowl_mdns, 'Bowling', ['C2:C' end_cell]);
    xlswrite(excel_out, bowl_runs, 'Bowling', ['D2:D' end_cell]);
    xlswrite(excel_out, bowl_wkts, 'Bowling', ['E2:E' end_cell]);
    xlswrite(excel_out, bowl_averages, 2, ['F2:F' end_cell]);
    xlswrite(excel_out, bowl_srs, 'Bowling', ['G2:G' end_cell]);
    xlswrite(excel_out, bowl_ers, 'Bowling', ['H2:H' end_cell]);
    xlswrite(excel_out, bowl_bfs, 'Bowling', ['I2:I' end_cell]);
end