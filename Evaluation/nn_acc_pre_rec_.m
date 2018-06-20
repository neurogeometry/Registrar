
clc
clear

% TO FIND THE K NEAREST NEIGHBOUR ::

% FORMULA :: 
        % Class = knnclassify(Test, Training, Group, k, distance, rule)

% INPUTS :-
        % TEST SET, TRAINING SET, GROUP, VALUE OF K, DISTANCE, RULE
        
% TEST, TRAINING - MATRIX
% GROUP - GROUPING OF ROWS
% K - NUMBER OF NEAREST NEIGHBOR USED IN CLASSIFICATION
% DISTANCE - EUCLIDEAN
% RULE - NEAREST
        
%OUTPUT :-
        % CLASS MATRIX WHICH SHOWS THE NEAREST NEIGHBORS OF EACH ROW ie.,
        % ROW 1 OF TEST IS CLOSEST TO ROW 3 OF TRAINING
        
load test.txt;
load train.txt;

% load test.txt, train.txt --> WE CAN USE TEST AND TRAINING SET DATA HERE

% TO FIND THE ROW AND COLUMN TEST AND TRAINING SET DATA :: 

[test_row test_col] = size(test);
[train_row train_col] = size(train);


for i = 1 : 300
     new_test(i, :) = test(i, :);
end 

[new_test_row new_test_col] = size(new_test);

group = [ 1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120;121;122;123;124;125;126;127;128;129;130;131;132;133;134;135;136;137;138;139;140;141;142;143;144;145;146;147;148;149;150;151;152;153;154;155;156;157;158;159;160;161;162;163;164;165;166;167;168;169;170;171;172;173;174;175;176;177;178;179;180;181;182;183;184;185;186;187;188;189;190;191;192;193;194;195;196;197;198;199;200;201;202;203;204;205;206;207;208;209;210;211;212;213;214;215;216;217;218;219;220;221;222;223;224;225;226;227;228;229;230;231;232;233;234;235;236;237;238;239;240;241;242;243;244;245;246;247;248;249;250;251;252;253;254;255;256;257;258;259;260;261;262;263;264;265;266;267;268;269;270;271;272;273;274;275;276;277;278;279;280;281;282;283;284;285;286;287;288;289;290;291;292;293;294;295;296;297;298;299;300 ];
     
disp(' RESULT ');
disp(' K  NEAREST NEIGHBOR :: ');

class = knnclassify( new_test, train, group, 1, 'euclidean', 'nearest' );

for i = 1 : 300
    fprintf( ' \n ROW %d OF TESTING SET IS CLOSEST TO ROW %d OF THE TRAINING SET ', i, class(i) );
end

fprintf( '\n \n' );

fprintf( '\n \n' );

true_positive = 0;
true_negative = 0;
false_positive = 0;
false_negative = 0;

for i = 1 : 300
    if( ( train(class(i), 11) == 4 ) && ( test(i, 11) == 4 ) )
        true_positive = true_positive + 1;
        elseif( ( train(class(i), 11) == 2 ) && ( test(i, 11) == 2 ) )
            true_negative = true_negative + 1;
            elseif( ( train(class(i), 11) == 2 ) && ( test(i, 11) == 4 ) )
                false_positive = false_positive + 1;
                else
                    false_negative = false_negative + 1;
    end    
end

positive = true_positive + true_negative;
negative = false_positive + false_negative;
true_positive;
true_negative;
false_positive;
false_negative;

accuracy = ( true_positive + true_negative ) / ( true_positive + true_negative + false_positive + false_negative )
precision = true_positive / ( true_positive + false_positive )
recall = true_positive / ( true_positive + false_negative )

mat = [ true_positive true_negative false_positive false_negative accuracy precision recall ]

plot( mat )

fprintf( '\n \n' );