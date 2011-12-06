PKG_NAME=$(shell grep Package ../DESCRIPTION | sed -e 's/.*: //')

index.html: index.html.in configure ../DESCRIPTION
	./configure

configure: configure.in
	autoconf

include $(OMEGA_HOME)/R/Config/RPackageWebConfig

WEB_DIR=www.omegahat.org:/home3/WebSites/Omega/$(PKG_NAME)

shipIndex: aux
	$(MAKE) -C ../inst/doc ship

aux: index.html ../demo/synthetic.R
	scp index.html $(WEB_DIR)/index.html
	scp ../demo/synthetic.R $(WEB_DIR)/demo
	-scp ../Changes ../FAQ.html $(WEB_DIR)

indexInstall indexOnly: index.html
	scp index.html $(WEB_DIR)/index.html
