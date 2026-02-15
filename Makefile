native: check-env
	cargo make build-versioned

check-env:
	@command -v cargo-make >/dev/null 2>&1 || { \
		echo "\033[0;31mError: 'cargo-make' is not installed.\033[0m"; \
		echo "\033[0;32mPlease install it by running\033[0m:\n\033[1m$$ cargo install cargo-make\033[0m"; \
		exit 1; \
	}

windows: check-env
	cargo make build-windows

wasm: check-env
	cargo make build-wasm
