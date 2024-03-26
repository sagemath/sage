# Prints gap kernel info without initializing gap.
#
# Usage:
#
#   $ gap -r --systemfile gap_kernel_info.g
#   KERNEL_VERSION=4.12.2
#   GAP_ARCHITECTURE=x86_64-unknown-linux-gnu-default64-kv8
#   GAP_ROOT_PATHS=/usr/lib64/gap/;/usr/share/gap/

KernelInfo := KERNEL_INFO();

PRINT_TO( "*stdout*", "KERNEL_VERSION=", KernelInfo.KERNEL_VERSION, "\n");
PRINT_TO( "*stdout*", "GAP_ARCHITECTURE=", KernelInfo.GAP_ARCHITECTURE, "\n");

p := KernelInfo.GAP_ROOT_PATHS;
PRINT_TO( "*stdout*", "GAP_ROOT_PATHS=", p[1] );
for i in [2 .. LENGTH(p)] do
  PRINT_TO( "*stdout*", ";", p[i]);
od;
PRINT_TO( "*stdout*", "\n");

QuitGap();
