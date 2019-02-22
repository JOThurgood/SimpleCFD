module welcome

  implicit none

  public

  contains

  subroutine welcome_msg
  
    call execute_command_line("clear")

    print *, '++++++++++++ooooooooooooooosssyyyyyyyyyyyyyyysssssssssyyyyyyyyyyyyyyyy'
    print *, 'ooo++++++++++++++/+oooooooossssyyyyyyyyyyyyyyyyssssssyysssssyyyyyyyyyy'
    print *, 'oooooo+++++/:--....-/ooosossssyyyyyyyyyyyyyyyyyyssyyso/-::/+sssssyyhhy'
    print *, 'ssssooooo/:.-:///:-..:+sosssssyyyyyyyyyyyyyyyyyyyss/-...-----:::oyyhhy'
    print *, 'ysssssss+--/osyys+/---:ssssoooo+++/+++//++osysysso:...-:+oso/:.-oyyhhy'
    print *, 'yyyyyyso+/+shdddhyo//++yyyss/-.```````````..-/osyy+/:-/oyhhhyo:-:syyyy'
    print *, 'yyyyyys+osshdmmddhyosyyyyyss+-.````````````-/oyyhyyoosyyyyddhyo::+syyy'
    print *, 'hyyyyyssyhhdmNNdhhhhhhhhhys+:.````````````.-/oyyhhhyhhhyssddhys+/+osss'
    print *, 'hyyyyysyhdmmmNmdddhhhhhhhys/-``````````````.:sssyyyhhhhyoshmdhhsooosss'
    print *, 'yyhyyysyddmNmmhddhhhhhhhhyy+:``````````````.-+ssssssyyhhyyyhddddhsoooo'
    print *, 'hhhyyys+yddmdhyhhhyhhdhhhyyo-```````.``````.:oyyyyysssyhhyyssyhdds++oo'
    print *, 'hhhyyss//syssossyyhhhdhhhys/:``````...`.```-oyyhyyysssosyysoo+oyy+/+++'
    print *, 'yyyyyo:---:++//osyhhhdhhhhs:-.``````....```:oyhhyyyyyysoso+//://:---:/'
    print *, 'yyyyys/-.------:oyhhhddddhy+-....````````../yhhhhhyhyyss/:----......-:'
    print *, 'yyyyso:-.......:oyhdddddddho-.`.`````````../sdhdhhyyyyo/-....`.......:'
    print *, 'yyyys+-.......--:ohdddmdddho:`````````````.:ohdddhhhys/-.``.........-:'
    print *, 'syyso/-....-----:+yhddmmmdho-.`````````````:ydddddhys/-.````........-:'
    print *, 'sssso:-...-------:ohddmmmdd+..`````````````/hddddddhy:...`.......`...-'
    print *, 'sss+/:...--------:+hdddddmh+...`.``````````:syhhhdddy-...............-'
    print *, 'so+:-.------------/hdddmmmds...``````.````.+yydddddds-...............-'
    print *, 'o+:-..--------:::--sddmddmmh...-``....````-ydddhdmdh+.....---........-'
    print *, 'o+:-------::://:---smdmmNmdh-..-.......```-yhdNNmdmd+....-----........'
    print *, 'o+:------:::://::--+mmddmdhd+-.........`.`:yhyhddddh+....-----.......-'
    print *, 'o/:----::::-:://:-.:yddhhhhho-......````../yyssyyyy/-...-::------.....'
    print *, 'o/:--::/+o+:::::::--:oyhhhhys-.....``````.osyhyhyo:...--:::::::::-....'
    print *, '+/:--:/+syys//::::--.-+shhhyo-........`..-osyhhh+-...-:////:::///:-...'
    print *, '/----:/+syhhys+++/:---:shhhy+-............+syhhs:..---/osssooo+/:--.--'
    print *, ':--::/+osyhddhyyys+:--:ohhhy+-........``../syys+--::/+syhhhhyys+/:---:'
    print *, ':-://+ossyhdmdhhhhys+/:oyso+:-....---..`.--/osyo:/+shhhhhhhhhyso+/::::'
    print *, '--::/+osyhddmmmdddhhyo++/:---...--:::-......-/+o/oshhhhhddddhyyssso+//'
    print *, ':://++osyhdmmmmmmdddddy+--..-:::++++++/:-.....-/syhhddmmmmmmddddhyso++'
    print *, '://++osyhddmmmmNmmmmmmms:--/ossyyyssssso++:...-oddmmNNNNNmmmmmddhyysoo'
    print *, '/++osyhhdddmmNNNNNNNNNNd+:+hyyyyyyyyyyyyyyy/-:+hmNNNNNNNmmmmmmddhso++o'
    print *, 'ossyhhddmmmmNNNNNNNNNNmmy/shhhyhhdddhyhddhho:ohNmmNNNNNNNNNmmddhyo+//+'
    print *, 'yyhhddmmmmNNNNNNNNNNNNNNNysdmmmmmdddmNNmmmdoodmmNmNNmNNNNNNmmddhyso+/+'
    print *, 'yyhdmmmNNNNNNNNNNNNMNMNNdmdddmmmmdddmmNNmdyhdmmmmNNNNNNNNNNNmmdhhyso+/'
    print *, 'yhhddmmNNNNNNNMMMMMMMNNNNNNmNmmNNmmmmmmNNmdmmmNNNNNNNMMNNNNNmmddhysso+'
    print *, 'yyhhdmmNNNNNMMMMMMMMMMMMNNNMNNNNNNNNNNNNNNNNNNNNNNNNMMMNNNNNmmdhhysso+'
    print *, 'yyhddmmNNNNNMMMMMMMMMMMMNMMMNNNNNNNNNNNNNNNNNNMNNMMMMMMNNNNNmddhhyso++'
    print *, 'syddmmmNNNNNMMMMMMMMMMMMMMMMMNNNNNNNNNNNNNNNNNNMMMMMMMMNNNNmmddddhysoo'
    print *, 'hddmmmNNNNNMMMMMMMMMMMMMMMMMMNNNNNNNNNNNNNNNNNNMMMMMMMMMNNNmmmmmmdhyys'
    print *, 'dmmmmNNNNMMMMMMMNNNNNNMMMNNNNNNNNNNNNNNNNNNNNNNMMMMMMMMNNNNNNmmmmmdhys'
    print *, 'mmmmmNNNNNMMMMMNNNNNNNNNMNNNNNNNNNNNNNNNNNNNNNNNNMMMMMMMNNNNNNmmmmmdhy'
    print *, 'mmmmNNNNNMMMMMMNNNNNNNNNNNNNNNNNNNNNNNNmmmNNNNNNNNNMMMMMMNNNNNNNmmmddh'
    print *, 'mmNNNNNNMMMMMMMNNNNNNNNNNNNNNNNNNNNNNNNmmmmmmNNNNNNNMMMMMNNNNNNNmmmddh'
    print *, 'mNNNNMNMMMMMMMNNNNNNNNNNNNNNNNNNNNNNNNNNNNNmmmNNNNNNNNMMMNNNNNNNmmmddd'
    print *, 'NNNMMMMNMMMMMMNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNmmNNNNNNNMMMMNNNNNNmmmmddd'
    print *, 'NNNNNNNNNNMMMNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMMMMMMNNNNNmmmmddd'
    print *, 'NNNNMMNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMMMMMMMMNNNmmmmmmd'
    print *, 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMMMMMMMMMNNNNmmmmmm'
    print *, 'NNNNNNNNNMNNNNNNNNNNNNNNNNNNNNNmmNNNNNNNNNNNNNNNNMMMMMMMMMMMMNNNNmmmmm'

    print *,''
    print *,''
    print *,''

    print *,'Welcome to CCAPS (Cell-centered approximate projection solver)'
    print *,''
    print *,'Using variant CCAPS_atm_const (subtracted 1d hydrostatic background)'
    print *,'******************************************************************'
    print *,'Launching code'
  end subroutine welcome_msg

end module welcome
