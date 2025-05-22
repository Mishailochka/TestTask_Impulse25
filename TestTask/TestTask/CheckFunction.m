%% Чистка
clc; clear;
close all;

%% Чтение из файла
% Чтение в ASCII-массив
    fid = fopen('BER_from_Cpp.txt', 'r');
    i = 1;
    stop = 0;
    while stop ~= 1
        element = fscanf(fid, '%c', 1); % Чтение одного символа
        
        if ~isempty(element)
            NumText(i) = double(element);
        else
            disp("Файл прочитан")
            stop = 1;
        end
        i = i+1;
    end
    fclose(fid);

% Разделение текста на строчки
    EndOfStr = find(NumText == 10);
    NumTextVAR = NumText(1:EndOfStr(1)-2);
    NumTextBER = NumText(EndOfStr(1)+1:EndOfStr(2)-2);
    NumTextConst = NumText(EndOfStr(2)+1:EndOfStr(3)-2);

% Обработка VAR
    EndOfCol = find(NumTextVAR == 9);
    Ncol = length(EndOfCol);
    VAR = NaN(Ncol, 1);

    EdgesCol = [0 EndOfCol];
    for i = 1 : Ncol
        str = NumTextVAR(EdgesCol(i)+1 : EdgesCol(i+1)-1);
        VAR(i) = str2double(char(str));
    end

% Обработка BER
    EndOfCol = find(NumTextBER == 9);
    Ncol = length(EndOfCol);
    BER = NaN(Ncol, 1);

    EdgesCol = [0 EndOfCol];
    for i = 1 : Ncol
        str = NumTextBER(EdgesCol(i)+1 : EdgesCol(i+1)-1);
        BER(i) = str2double(char(str));
    end

% Обработка Constellation
    EndOfCol = find(NumTextConst == 41);
    CommaPlace = find(NumTextConst == 44);
    Ncol = length(EndOfCol);
    Const = NaN(Ncol, 1);

    EdgesCol = [-1 EndOfCol];
    for i = 1 : Ncol
        str1 = NumTextConst(EdgesCol(i)+3 : CommaPlace(i)-1);
        str2 = NumTextConst(CommaPlace(i)+1 : EdgesCol(i+1)-1);
        Const(i) = complex( str2double(char(str1)), str2double(char(str2)) );
    end

%% Сравниваем с теоретической
% berawgn() считает кривую для SNR, поэтому надо для начала VAR перевечти в
% SNR. Для этого как раз-таки и нужно сигнальное созвездие, которое мы
% передавали вместе с BER и VAR из C++

    Es = mean(abs(Const).^2);
    Eb = Es / log2(length(Const));
    
    N0 = 2*VAR;
    SNR = 10*log10(Eb./N0);
    
    BERth = berawgn(SNR, "qam", length(Const));
    
    figure(1)
        semilogy(VAR, BER, 'LineWidth', 2, 'Color', 'blue', 'LineStyle', '-', 'DisplayName', 'C++ modeling');
        hold on;
        semilogy(VAR, BERth, 'LineWidth', 2, 'Color', 'red', 'LineStyle', '--', 'DisplayName', 'MatLab theory');
    
        grid on; grid minor;
        legend('Location', 'southeast');

        xlabel("VAR")
        ylabel("BER")