#!/bin/bash
# ==============================================================================
# GPU UTILITIES
# ==============================================================================
# GPU detection, configuration, and setup for pipeline acceleration
# Sourced by: modules_loader.sh, gpu_prep.sh
# ==============================================================================

# Guard against double-sourcing
[[ "${GPU_UTILS_SOURCED:-}" == "true" ]] && return 0
export GPU_UTILS_SOURCED="true"

# ==============================================================================
# GPU CONFIGURATION
# ==============================================================================

GPU_AVAILABLE="false"
GPU_COUNT=0
GPU_MEMORY_MB=0
CUDA_VERSION=""
CUDA_READY="false"
GPU_VENDOR=""
DISTRO=""
DISTRO_VERSION=""

# ==============================================================================
# LOGGING CONFIGURATION
# ==============================================================================

GPU_LOG_DIR="${GPU_LOG_DIR:-${SCRIPT_DIR:-$(pwd)}/logs}"
GPU_LOG_FILE="${GPU_LOG_FILE:-$GPU_LOG_DIR/gpu_prep_$(date +%Y%m%d_%H%M%S).log}"
GPU_LOG_ENABLED="${GPU_LOG_ENABLED:-true}"
GPU_LOG_LEVEL="${GPU_LOG_LEVEL:-INFO}"  # DEBUG, INFO, WARN, ERROR

# Ensure log directory exists
mkdir -p "$GPU_LOG_DIR" 2>/dev/null || true

# ==============================================================================
# LOGGING FUNCTIONS
# ==============================================================================

# Get timestamp for logging
_get_timestamp() {
	date "+%Y-%m-%d %H:%M:%S"
}

# Write to log file
_write_log() {
	local level="$1"
	local message="$2"
	if [[ "$GPU_LOG_ENABLED" == "true" && -n "$GPU_LOG_FILE" ]]; then
		echo "[$(_get_timestamp)] [$level] $message" >> "$GPU_LOG_FILE" 2>/dev/null || true
	fi
}

# Log level check (returns 0 if should log)
_should_log() {
	local level="$1"
	case "$GPU_LOG_LEVEL" in
		DEBUG) return 0 ;;
		INFO)  [[ "$level" != "DEBUG" ]] && return 0 ;;
		WARN)  [[ "$level" == "WARN" || "$level" == "ERROR" ]] && return 0 ;;
		ERROR) [[ "$level" == "ERROR" ]] && return 0 ;;
	esac
	return 1
}

# Logging stubs (if not already defined by another module)
if ! declare -f log_info &>/dev/null; then
	log_debug() {
		if _should_log "DEBUG"; then
			echo "[DEBUG] $*"
			_write_log "DEBUG" "$*"
		fi
	}
	log_info() {
		if _should_log "INFO"; then
			echo "[INFO] $*"
			_write_log "INFO" "$*"
		fi
	}
	log_warn() {
		if _should_log "WARN"; then
			echo "[WARN] $*" >&2
			_write_log "WARN" "$*"
		fi
	}
	log_error() {
		if _should_log "ERROR"; then
			echo "[ERROR] $*" >&2
			_write_log "ERROR" "$*"
		fi
	}
	log_step() {
		echo "=== $* ==="
		_write_log "STEP" "$*"
	}
fi

# Log system info at start
_log_system_info() {
	_write_log "INFO" "========================================"
	_write_log "INFO" "GPU Utils initialized"
	_write_log "INFO" "Host: $(hostname)"
	_write_log "INFO" "User: $(whoami)"
	_write_log "INFO" "Kernel: $(uname -r)"
	_write_log "INFO" "WSL: $(is_wsl && echo 'yes' || echo 'no')"
	_write_log "INFO" "========================================"
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

check_root() {
	if [[ $EUID -eq 0 ]]; then
		log_error "This script should not be run as root"
		exit 1
	fi
}

detect_distro() {
	if [[ -f /etc/os-release ]]; then
		. /etc/os-release
		DISTRO="$ID"
		DISTRO_VERSION="$VERSION_ID"
	else
		DISTRO="unknown"
		DISTRO_VERSION="unknown"
	fi
	export DISTRO DISTRO_VERSION
}

detect_gpu_type() {
	GPU_VENDOR=""
	# Check via lspci first (native Linux)
	if lspci 2>/dev/null | grep -qi "nvidia"; then
		GPU_VENDOR="nvidia"
	elif lspci 2>/dev/null | grep -qi "amd.*radeon\|amd.*vega\|amd.*navi"; then
		GPU_VENDOR="amd"
	elif lspci 2>/dev/null | grep -qi "intel.*graphics\|intel.*iris"; then
		GPU_VENDOR="intel"
	# WSL2 fallback: check for nvidia-smi directly (GPU not visible via lspci in WSL2)
	elif command -v nvidia-smi &>/dev/null && nvidia-smi &>/dev/null; then
		GPU_VENDOR="nvidia"
		log_info "WSL2 detected - using nvidia-smi for GPU detection"
	fi
	export GPU_VENDOR
	log_info "Detected GPU vendor: ${GPU_VENDOR:-none}"
}

# Check if running in WSL
is_wsl() {
	[[ -f /proc/version ]] && grep -qi "microsoft\|wsl" /proc/version 2>/dev/null
}

install_base_packages() {
	log_step "Installing base packages"
	
	# Skip kernel headers on WSL2 (not available/needed)
	local kernel_pkg=""
	if ! is_wsl; then
		kernel_pkg="dkms linux-headers-$(uname -r)"
	else
		log_info "WSL2 detected - skipping kernel headers (not needed)"
	fi
	
	case "$DISTRO" in
		ubuntu|debian)
			sudo apt-get update -qq
			sudo apt-get install -y wget curl software-properties-common gnupg \
				build-essential $kernel_pkg \
				htop sysstat iotop
			;;
		fedora|rhel|centos|rocky|almalinux)
			if ! is_wsl; then
				kernel_pkg="kernel-devel kernel-headers dkms"
			fi
			sudo dnf install -y wget curl gnupg2 \
				gcc make $kernel_pkg \
				htop sysstat iotop
			;;
		*)
			log_warn "Unknown distro: $DISTRO - skipping base packages"
			;;
	esac
}

# ==============================================================================
# GPU DETECTION
# ==============================================================================

detect_gpu() {
	GPU_AVAILABLE="false"
	GPU_COUNT=0
	CUDA_READY="false"
	
	# Check for NVIDIA GPU via nvidia-smi
	if command -v nvidia-smi &>/dev/null; then
		local gpu_info
		gpu_info=$(nvidia-smi --query-gpu=count,memory.total --format=csv,noheader,nounits 2>/dev/null)
		if [[ -n "$gpu_info" ]]; then
			GPU_AVAILABLE="true"
			GPU_COUNT=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | wc -l)
			GPU_MEMORY_MB=$(echo "$gpu_info" | head -1 | cut -d',' -f2 | tr -d ' ')
			CUDA_VERSION=$(nvidia-smi 2>/dev/null | grep -oP 'CUDA Version: \K[0-9.]+' | head -1)
			
			# Check if CUDA toolkit is ready
			if command -v nvcc &>/dev/null; then
				CUDA_READY="true"
			fi
		fi
	fi
	
	export GPU_AVAILABLE GPU_COUNT GPU_MEMORY_MB CUDA_VERSION CUDA_READY
}

# Check if GPU is available
has_gpu() {
	[[ "$GPU_AVAILABLE" == "true" ]]
}

# Check if CUDA is ready for use
is_cuda_ready() {
	[[ "$CUDA_READY" == "true" ]]
}

# Log GPU status
log_gpu_status() {
	if has_gpu; then
		log_info "[GPU] Detected $GPU_COUNT GPU(s), ${GPU_MEMORY_MB}MB memory, CUDA: $CUDA_VERSION"
		if is_cuda_ready; then
			log_info "[GPU] CUDA toolkit ready (nvcc available)"
		else
			log_warn "[GPU] CUDA toolkit not installed - run setup_cuda for GPU acceleration"
		fi
	else
		log_info "[GPU] No GPU detected - using CPU only"
	fi
}

# ==============================================================================
# CUDA SETUP AND PREPARATION
# ==============================================================================

# Check CUDA prerequisites
check_cuda_prerequisites() {
	local missing=()
	
	# Check for nvidia driver
	if ! command -v nvidia-smi &>/dev/null; then
		missing+=("nvidia-driver")
	fi
	
	# Check for nvcc (CUDA compiler)
	if ! command -v nvcc &>/dev/null; then
		missing+=("cuda-toolkit")
	fi
	
	# Check for cuDNN headers
	if [[ ! -f "/usr/include/cudnn.h" && ! -f "/usr/local/cuda/include/cudnn.h" ]]; then
		missing+=("cudnn")
	fi
	
	if [[ ${#missing[@]} -gt 0 ]]; then
		echo "${missing[*]}"
		return 1
	fi
	return 0
}

# Setup CUDA environment
setup_cuda() {
	log_step "Setting up CUDA environment"
	
	# Detect GPU first
	detect_gpu
	
	if ! has_gpu; then
		log_error "No NVIDIA GPU detected. Cannot setup CUDA."
		return 1
	fi
	
	# Check what's missing
	local missing
	missing=$(check_cuda_prerequisites)
	
	if [[ -z "$missing" ]]; then
		log_info "[CUDA] All prerequisites already installed"
		configure_cuda_env
		return 0
	fi
	
	log_info "[CUDA] Missing components: $missing"
	log_info "[CUDA] Attempting installation..."
	
	# Detect package manager
	local pkg_mgr=""
	if command -v apt-get &>/dev/null; then
		pkg_mgr="apt"
	elif command -v yum &>/dev/null; then
		pkg_mgr="yum"
	elif command -v dnf &>/dev/null; then
		pkg_mgr="dnf"
	elif command -v conda &>/dev/null; then
		pkg_mgr="conda"
	fi
	
	case "$pkg_mgr" in
		apt)
			install_cuda_apt
			;;
		yum|dnf)
			install_cuda_yum
			;;
		conda)
			install_cuda_conda
			;;
		*)
			log_error "Unsupported package manager. Please install CUDA manually."
			log_info "See: https://developer.nvidia.com/cuda-downloads"
			return 1
			;;
	esac
	
	configure_cuda_env
	detect_gpu  # Re-detect after installation
	log_gpu_status
}

# Install CUDA via apt (Ubuntu/Debian)
install_cuda_apt() {
	log_info "[CUDA] Installing via apt..."
	
	# Add NVIDIA package repository
	if [[ ! -f "/etc/apt/sources.list.d/cuda.list" ]]; then
		log_info "[CUDA] Adding NVIDIA repository..."
		wget -qO - https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/3bf863cc.pub | \
			sudo gpg --dearmor -o /usr/share/keyrings/cuda-archive-keyring.gpg 2>/dev/null || true
	fi
	
	sudo apt-get update -qq
	sudo apt-get install -y cuda-toolkit nvidia-cuda-toolkit 2>/dev/null || \
		log_warn "apt install failed - try manual CUDA installation"
}

# Install CUDA via yum/dnf (RHEL/CentOS/Fedora)
install_cuda_yum() {
	log_info "[CUDA] Installing via yum/dnf..."
	sudo yum install -y cuda-toolkit 2>/dev/null || \
	sudo dnf install -y cuda-toolkit 2>/dev/null || \
		log_warn "yum/dnf install failed - try manual CUDA installation"
}

# Install CUDA via conda
install_cuda_conda() {
	log_info "[CUDA] Installing via conda..."
	conda install -y -c nvidia cuda-toolkit cudnn 2>/dev/null || \
		log_warn "conda install failed - try: conda install -c nvidia cuda-toolkit"
}

# Configure CUDA environment variables
configure_cuda_env() {
	log_info "[CUDA] Configuring environment..."
	
	# Find CUDA installation
	local cuda_home=""
	for path in /usr/local/cuda /usr/local/cuda-* /opt/cuda; do
		if [[ -d "$path" ]]; then
			cuda_home="$path"
			break
		fi
	done
	
	if [[ -n "$cuda_home" ]]; then
		export CUDA_HOME="$cuda_home"
		export PATH="$CUDA_HOME/bin:$PATH"
		export LD_LIBRARY_PATH="$CUDA_HOME/lib64:${LD_LIBRARY_PATH:-}"
		log_info "[CUDA] CUDA_HOME=$CUDA_HOME"
	fi
	
	# cuDNN paths
	if [[ -d "/usr/local/cudnn" ]]; then
		export LD_LIBRARY_PATH="/usr/local/cudnn/lib64:$LD_LIBRARY_PATH"
	fi
}

# Configure CUDA environment and add to bashrc (standalone version)
configure_cuda_env_standalone() {
	log_info "Configuring CUDA environment..."
	
	local cuda_home=""
	for path in /usr/local/cuda /usr/local/cuda-* /opt/cuda; do
		if [[ -d "$path" ]]; then
			cuda_home="$path"
			break
		fi
	done
	
	if [[ -n "$cuda_home" ]]; then
		# Add to bashrc if not already present
		if ! grep -q "CUDA_HOME" ~/.bashrc 2>/dev/null; then
			cat >> ~/.bashrc << EOF

# CUDA Configuration
export CUDA_HOME="$cuda_home"
export PATH="\$CUDA_HOME/bin:\$PATH"
export LD_LIBRARY_PATH="\$CUDA_HOME/lib64:\$LD_LIBRARY_PATH"
EOF
			log_info "Added CUDA to ~/.bashrc"
		fi
		
		export CUDA_HOME="$cuda_home"
		export PATH="$CUDA_HOME/bin:$PATH"
		export LD_LIBRARY_PATH="$CUDA_HOME/lib64:${LD_LIBRARY_PATH:-}"
	fi
}

# ==============================================================================
# NVIDIA GPU SETUP
# ==============================================================================

install_nvidia_driver() {
	log_step "Installing NVIDIA driver v${NVIDIA_DRIVER_VERSION:-535}"
	
	case "$DISTRO" in
		ubuntu|debian)
			# Add NVIDIA repository
			local ubuntu_ver=$(lsb_release -rs 2>/dev/null | sed 's/\.//' || echo "2204")
			local keyring_url="https://developer.download.nvidia.com/compute/cuda/repos/ubuntu${ubuntu_ver}/x86_64/cuda-keyring_1.1-1_all.deb"
			
			if [[ ! -f /usr/share/keyrings/cuda-archive-keyring.gpg ]]; then
				log_info "Adding NVIDIA CUDA repository..."
				wget -q "$keyring_url" -O /tmp/cuda-keyring.deb 2>/dev/null || \
					wget -q "https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb" -O /tmp/cuda-keyring.deb
				sudo dpkg -i /tmp/cuda-keyring.deb
				rm -f /tmp/cuda-keyring.deb
				sudo apt-get update -qq
			fi
			
			sudo apt-get install -y "nvidia-driver-${NVIDIA_DRIVER_VERSION:-535}" \
				"nvidia-utils-${NVIDIA_DRIVER_VERSION:-535}" 2>/dev/null || \
			sudo apt-get install -y nvidia-driver nvidia-utils
			;;
		fedora|rhel|centos|rocky|almalinux)
			sudo dnf config-manager --add-repo \
				https://developer.download.nvidia.com/compute/cuda/repos/rhel9/x86_64/cuda-rhel9.repo 2>/dev/null || true
			sudo dnf install -y nvidia-driver nvidia-settings
			;;
	esac
}

install_cuda_toolkit() {
	log_step "Installing CUDA Toolkit"
	
	case "$DISTRO" in
		ubuntu|debian)
			sudo apt-get install -y cuda-toolkit nvidia-cuda-toolkit 2>/dev/null || \
				sudo apt-get install -y cuda
			;;
		fedora|rhel|centos|rocky|almalinux)
			sudo dnf install -y cuda-toolkit
			;;
	esac
	
	# Configure environment
	configure_cuda_env 2>/dev/null || configure_cuda_env_standalone
}

install_cudnn() {
	log_step "Installing cuDNN"
	
	case "$DISTRO" in
		ubuntu|debian)
			sudo apt-get install -y libcudnn8 libcudnn8-dev 2>/dev/null || \
				log_warn "cuDNN not available via apt - install manually from NVIDIA"
			;;
		fedora|rhel|centos|rocky|almalinux)
			sudo dnf install -y cudnn 2>/dev/null || \
				log_warn "cuDNN not available via dnf - install manually from NVIDIA"
			;;
	esac
}

install_nvidia_monitoring() {
	log_step "Installing NVIDIA monitoring tools"
	
	pip3 install --user --quiet nvidia-ml-py gpustat 2>/dev/null || \
		pip install --user --quiet nvidia-ml-py gpustat 2>/dev/null || \
		log_warn "Could not install Python GPU monitoring packages"
}

setup_nvidia() {
	log_step "Setting up NVIDIA GPU"
	install_nvidia_driver
	install_cuda_toolkit
	install_cudnn
	install_nvidia_monitoring
}

# ==============================================================================
# AMD GPU SETUP
# ==============================================================================

setup_amd() {
	log_step "Setting up AMD GPU (ROCm)"
	
	case "$DISTRO" in
		ubuntu|debian)
			# Add ROCm repository
			wget -q -O - https://repo.radeon.com/rocm/rocm.gpg.key | sudo apt-key add - 2>/dev/null || true
			echo "deb [arch=amd64] https://repo.radeon.com/rocm/apt/debian/ ubuntu main" | \
				sudo tee /etc/apt/sources.list.d/rocm.list
			sudo apt-get update -qq
			sudo apt-get install -y rocm-smi-lib rocm-dev 2>/dev/null || \
				log_warn "ROCm installation failed - check AMD support for your GPU"
			;;
		fedora|rhel|centos|rocky|almalinux)
			sudo dnf install -y rocm-smi rocm-dev 2>/dev/null || \
				log_warn "ROCm installation failed"
			;;
	esac
}

# ==============================================================================
# INTEL GPU SETUP
# ==============================================================================

setup_intel() {
	log_step "Setting up Intel GPU"
	
	case "$DISTRO" in
		ubuntu|debian)
			sudo apt-get install -y intel-gpu-tools intel-media-va-driver \
				vainfo 2>/dev/null || log_warn "Intel GPU tools installation failed"
			;;
		fedora|rhel|centos|rocky|almalinux)
			sudo dnf install -y intel-gpu-tools libva-intel-driver 2>/dev/null || \
				log_warn "Intel GPU tools installation failed"
			;;
	esac
}

# ==============================================================================
# VERIFICATION
# ==============================================================================

verify_installation() {
	log_step "Verifying GPU installation"
	
	local success=true
	
	case "$GPU_VENDOR" in
		nvidia)
			if command -v nvidia-smi &>/dev/null; then
				log_info "nvidia-smi: OK"
				nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
			else
				log_warn "nvidia-smi: NOT FOUND"
				success=false
			fi
			
			if command -v nvcc &>/dev/null; then
				log_info "CUDA compiler: OK ($(nvcc --version | grep release | awk '{print $5}' | tr -d ','))"
			else
				log_warn "CUDA compiler (nvcc): NOT FOUND - reboot may be required"
			fi
			;;
		amd)
			if command -v rocm-smi &>/dev/null; then
				log_info "rocm-smi: OK"
			else
				log_warn "rocm-smi: NOT FOUND"
				success=false
			fi
			;;
		intel)
			if command -v intel_gpu_top &>/dev/null; then
				log_info "intel_gpu_top: OK"
			else
				log_warn "intel_gpu_top: NOT FOUND"
				success=false
			fi
			;;
	esac
	
	# Check if reboot needed
	if [[ "$GPU_VENDOR" == "nvidia" ]]; then
		if ! lsmod | grep -q nvidia; then
			log_warn "NVIDIA kernel module not loaded - REBOOT REQUIRED"
			success=false
		fi
	fi
	
	$success
}

# ==============================================================================
# GPU-ACCELERATED TOOL DETECTION
# ==============================================================================

# Check if GPU-accelerated STAR is available
has_gpu_star() {
	has_gpu && command -v STARlong &>/dev/null
}

# Check for GPU-accelerated tools
check_gpu_tools() {
	log_step "Checking GPU-accelerated tools"
	
	local tools_found=0
	
	# RAPIDS cuML for ML acceleration
	if python -c "import cuml" 2>/dev/null; then
		log_info "[GPU-TOOL] RAPIDS cuML available"
		((tools_found++))
	fi
	
	# GPU-accelerated compression
	if command -v nvcomp &>/dev/null; then
		log_info "[GPU-TOOL] nvcomp (GPU compression) available"
		((tools_found++))
	fi
	
	# NVIDIA DALI for data loading
	if python -c "import nvidia.dali" 2>/dev/null; then
		log_info "[GPU-TOOL] NVIDIA DALI available"
		((tools_found++))
	fi
	
	if [[ $tools_found -eq 0 ]]; then
		log_info "[GPU-TOOL] No GPU-accelerated bioinformatics tools detected"
	fi
}

# ==============================================================================
# GPU-AWARE CONFIGURATION
# ==============================================================================

# Get optimal STAR genome load mode based on GPU/memory
get_star_genome_load() {
	if has_gpu && [[ $GPU_MEMORY_MB -gt 8000 ]]; then
		echo "LoadAndKeep"
	else
		echo "NoSharedMemory"
	fi
}

# Get optimal thread count considering GPU
get_optimal_threads() {
	local base_threads="${1:-$THREADS}"
	if has_gpu; then
		echo $((base_threads + 2))
	else
		echo "$base_threads"
	fi
}

# ==============================================================================
# MAIN FUNCTION FOR GPU PREPARATION
# ==============================================================================

gpu_prep_main() {
	log_step "GPU Preparation Script"
	log_info "Log file: $GPU_LOG_FILE"
	
	check_root
	detect_distro
	log_info "Detected: $DISTRO $DISTRO_VERSION"
	
	detect_gpu_type
	
	if [[ -z "$GPU_VENDOR" ]]; then
		log_warn "No dedicated GPU detected"
		log_info "Installing basic monitoring tools only..."
		install_base_packages
		return 0
	fi
	
	install_base_packages
	
	case "$GPU_VENDOR" in
		nvidia) setup_nvidia ;;
		amd)    setup_amd ;;
		intel)  setup_intel ;;
	esac
	
	echo ""
	verify_installation
	
	echo ""
	log_step "Setup Complete"
	log_info "To monitor GPU: nvidia-smi/rocm-smi"
	log_info "To use GPU in pipeline: source modules/logging/gpu_utils.sh"
	log_info "Log saved to: $GPU_LOG_FILE"
	
	if [[ "$GPU_VENDOR" == "nvidia" ]] && ! lsmod | grep -q nvidia; then
		echo ""
		log_warn ">>> REBOOT REQUIRED for GPU drivers to load <<<"
	fi
}

# ==============================================================================
# INITIALIZE GPU DETECTION
# ==============================================================================

# Initialize logging
_log_system_info

detect_gpu
configure_cuda_env 2>/dev/null || true

# Log GPU detection results
_write_log "INFO" "GPU_AVAILABLE=$GPU_AVAILABLE, GPU_COUNT=$GPU_COUNT, GPU_MEMORY_MB=$GPU_MEMORY_MB"
_write_log "INFO" "CUDA_VERSION=$CUDA_VERSION, CUDA_READY=$CUDA_READY"
