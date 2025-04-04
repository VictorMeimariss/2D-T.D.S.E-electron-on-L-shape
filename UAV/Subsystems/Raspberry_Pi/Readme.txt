This is the command you use in the pi to create the backup of the most significant files.

------> sudo tar -czvf pi_settings_backup.tar.gz   /etc/   /home/pi/   /boot/   /usr/local/   /var/lib/   /opt/

******if you want to Copy the /SDVNC folder as well use this command:
------> sudo tar -czvf pi_settings_backup_with_SDVNC.tar.gz   /etc/   /home/pi/   /home/SDVNC   /boot/   /usr/local/   /var/lib/   /opt/
******Having said that, due to this folder always changing, this should be happening once in a while not always, in order
for backup not taking that long

This command with your User_name and pi's ip is the one you use to copy the backup from the pi to your deskop automatically:

------> scp SDVNC@ip:~/pi_settings_backup.tar.gz "C:\Users\User_name\Desktop"

If any SD happens to burn:

To copy to pi

--->scp "C:\Users\User_name\Desktop\pi_settings_backup.tar.gz" new_pi_name@new_pi_ip:~/ 

This will overwrite system files with your backup.
----> sudo tar -xzvf pi_settings_backup.tar.gz -C /

then sudo reboot, if anything doesn't work well, search it up/enable/disable something







Useful commands:

To get files from the pi: scp SDVNC@ip:/home/SDVNC/path path/to/file example Desktop

To send files to the pi: scp path/to/file SDVNC@ip:/home/SDVNC/path

Path directly to pi "C:\Users\User_name\Documents\GitHub\UAV\Subsystems\Raspberry_Pi"

If you want to send folders use -r