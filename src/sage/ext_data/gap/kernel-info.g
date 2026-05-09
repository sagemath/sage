# Print the GAP version, architecture, and root paths on separate
# lines, without initializing the library. (That is, quickly.)
#
# Usage:
#
#   $ gap --systemfile kernel-info.g
#   4.15.1
#   riscv
#   /home/mjo/.gap/;/usr/lib/gap/;/usr/share/gap/
#
KernelInfo := KERNEL_INFO();
PRINT_TO( "*stdout*", KernelInfo.KERNEL_VERSION, "\n");
PRINT_TO( "*stdout*", KernelInfo.GAP_ARCHITECTURE, "\n");
p := KernelInfo.GAP_ROOT_PATHS;
PRINT_TO( "*stdout*", p[1] );
for i in [2 .. LENGTH(p)] do
  PRINT_TO( "*stdout*", ";", p[i]);
od;
PRINT_TO( "*stdout*", "\n");
QuitGap();
